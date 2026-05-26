# -*- coding: utf-8 -*-
"""
Background QgsTask for the habitat-free ASF Wildboar Connectivity plugin.

WORKFLOW (no core habitats needed):

    Outbreak point  ->  one Dijkstra from the origin disc
                    ->  cost-from-origin raster              (optional)
                    ->  continuous infection-risk raster     (exp(-cost/scale))
                    ->  one Circuit-theory solve             (origin -> AOI boundary)
                          via Circuitscape.jl if available
                          otherwise scipy.sparse fallback
                    ->  optional random-walk visit density

The AOI is determined automatically from the cost grid: cells with cost
above a multiple of the median become unreachable (effectively walls).
"""

import os
import tempfile

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.windows import (
    from_bounds,
    Window,
    transform as window_transform,
)
from rasterio.transform import xy as rio_xy

from qgis.core import (
    Qgis,
    QgsFeature,
    QgsField,
    QgsGeometry,
    QgsMessageLog,
    QgsPointXY,
    QgsProject,
    QgsRasterLayer,
    QgsTask,
    QgsVectorLayer,
)
from qgis.PyQt.QtCore import QVariant
from qgis.PyQt.QtGui import QColor
from qgis.PyQt.QtWidgets import QMessageBox

from .circuitscape_jl import (
    CircuitscapeJlError,
    is_julia_available,
    is_circuitscape_installed,
    run_one_to_all as run_circuitscape_jl,
)
from .layer_utils import ZOrder, add_wildboar_layer
from .resistance_editor import apply_modifications


LOG_TAG = "wildboar-asf"


class AsfConnectivityTask(QgsTask):
    """One-outbreak ASF connectivity analysis - no habitats input."""

    def __init__(self,
                 raster_path: str,
                 origin_xy: tuple,            # (x, y) in raster CRS
                 origin_radius_cells: int,    # disc radius in pixels (default 2)
                 fences: list,                # list of LineString geojson dicts
                 overpasses: list,            # list of (x, y) tuples
                 window_bounds: tuple,
                 crs_wkt: str,
                 options: dict):
        super().__init__("Wildboar ASF connectivity", QgsTask.CanCancel)
        self.raster_path = raster_path
        self.origin_xy = origin_xy
        self.origin_radius_cells = max(1, int(origin_radius_cells))
        self.fences = fences
        self.overpasses = overpasses
        self.window_bounds = window_bounds
        self.crs_wkt = crs_wkt
        self.options = options

        # Outputs
        self.cost_raster_path = None
        self.risk_raster_path = None
        self.current_raster_path = None
        self.walk_raster_path = None
        self.merged_lcp_polylines = []   # list[{"points", "traffic"}]
        self.forest_anchors = []         # list[((x, y), size)] - LCP targets
        self.n_lcps = 0                  # number of synthetic LCPs computed
        self.exception = None

        # Origin-snap tracking: if the user's click landed on an
        # impassable (NaN) cell, _origin_disc_mask snaps it to the
        # nearest passable cell and records the original + final
        # coordinates so the UI can show both.
        self.origin_was_snapped = False
        self.original_origin_xy = tuple(self.origin_xy)

        # Internal state
        self.nodata_mask = None
        self.aoi_mask = None
        self._cost_grid = None

    # =================================================================
    # Worker thread
    # =================================================================
    def run(self):
        try:
            arr, win_tf, crs, rows, cols = self._load_window()

            # Bake fences / overpasses into the working resistance grid.
            if self.fences or self.overpasses:
                arr = apply_modifications(
                    arr, win_tf, self.fences, self.overpasses,
                    fence_width_cells=2,
                    overpass_radius_cells=2,
                )
                QgsMessageLog.logMessage(
                    f"Applied {len(self.fences)} fence(s), "
                    f"{len(self.overpasses)} overpass(es) "
                    f"to the resistance window.",
                    LOG_TAG, Qgis.Info)

            # Cubic boar-resistance penalty: high-R pixels become near-walls.
            arr = self._apply_boar_resistance_penalty(arr)

            # ---- Rasterise the origin as a small disc on the grid ----
            # Pass the resistance array so an impassable click (highway,
            # lake, fence) snaps to the nearest passable cell instead of
            # failing.
            origin_cells = self._origin_disc_mask(arr.shape, win_tf, arr=arr)
            if origin_cells.size == 0:
                raise RuntimeError(
                    "No passable cells within ~2.5 km of the outbreak click. "
                    "Pick a point on land that is not entirely surrounded by "
                    "highways, large lakes, or fences.")
            QgsMessageLog.logMessage(
                f"Origin disc: {origin_cells.size} cell(s) at "
                f"({self.origin_xy[0]:.1f}, {self.origin_xy[1]:.1f}); "
                f"window {rows}x{cols}.",
                LOG_TAG, Qgis.Info)

            self.setProgress(10)
            if self.isCanceled():
                return False

            # ---- Single Dijkstra from origin -------------------------
            cost_grid, mcp = self._single_source_dijkstra(arr, origin_cells)
            self._cost_grid = cost_grid
            self.setProgress(25)
            if self.isCanceled():
                return False

            # ---- Auto AOI from the cost grid -------------------------
            self.aoi_mask = self._auto_aoi_mask(cost_grid, origin_cells,
                                                arr.shape)
            QgsMessageLog.logMessage(
                f"AOI: {int(self.aoi_mask.sum())} cells "
                f"({100 * self.aoi_mask.mean():.1f} % of window).",
                LOG_TAG, Qgis.Info)

            # ---- LCPs to AUTO-DETECTED FOREST destinations -----------
            if self.options.get("lcp", True):
                self._build_lcps_to_forests(
                    arr, win_tf, mcp,
                    n_max_targets=int(self.options.get("n_lcp_targets", 12)))
            self.setProgress(45)
            if self.isCanceled():
                return False

            # ---- Continuous infection-risk RASTER --------------------
            if self.options.get("risk", True):
                self.risk_raster_path = self._build_risk_raster(
                    cost_grid, win_tf)
            self.setProgress(50)
            if self.isCanceled():
                return False

            # ---- Pinchpoint raster (Circuitscape.jl preferred) -------
            if self.options.get("circuit", True):
                used_jl = False
                if self.options.get("use_circuitscape_jl", True) \
                        and is_julia_available() \
                        and is_circuitscape_installed():
                    try:
                        used_jl = self._try_circuitscape_jl(
                            arr, win_tf, origin_cells)
                        if used_jl:
                            QgsMessageLog.logMessage(
                                "Pinchpoints computed via Circuitscape.jl.",
                                LOG_TAG, Qgis.Success)
                    except CircuitscapeJlError as exc:
                        QgsMessageLog.logMessage(
                            f"Circuitscape.jl failed ({exc}); "
                            f"falling back to scipy solver.",
                            LOG_TAG, Qgis.Warning)
                        used_jl = False
                    except Exception as exc:
                        QgsMessageLog.logMessage(
                            f"Circuitscape.jl crashed ({exc}); "
                            f"falling back to scipy solver.",
                            LOG_TAG, Qgis.Warning)
                        used_jl = False
                if not used_jl:
                    self._single_source_current_scipy(
                        arr, win_tf, crs, origin_cells)
            self.setProgress(80)
            if self.isCanceled():
                return False

            self.setProgress(100)
            return True

        except Exception as exc:
            self.exception = exc
            QgsMessageLog.logMessage(
                f"Task failed: {exc}", LOG_TAG, Qgis.Critical)
            return False

    # =================================================================
    # Loading and pre-processing
    # =================================================================
    def _load_window(self):
        min_x, min_y, max_x, max_y = self.window_bounds
        with rasterio.open(self.raster_path) as src:
            crs = src.crs
            win = from_bounds(min_x, min_y, max_x, max_y, src.transform)
            # CRITICAL: snap the window to integer pixel boundaries so
            # the output rasters share the EXACT pixel grid with the
            # source resistance raster. Without this, the user's click
            # produces fractional col_off / row_off values (the click
            # rarely lands on a pixel corner), and every downstream
            # GeoTIFF ends up offset by some fraction of a pixel - the
            # ~17 m shift visible in QGIS when overlaying the pinchpoint
            # raster on the resistance raster.
            col_off = int(np.floor(win.col_off))
            row_off = int(np.floor(win.row_off))
            col_end = int(np.ceil(win.col_off + win.width))
            row_end = int(np.ceil(win.row_off + win.height))
            win = Window(col_off=col_off, row_off=row_off,
                         width=col_end - col_off,
                         height=row_end - row_off)
            win = win.intersection(Window(0, 0, src.width, src.height))
            if win.width <= 0 or win.height <= 0:
                raise RuntimeError(
                    "Analysis window does not overlap the resistance raster.")
            arr = src.read(1, window=win, boundless=True,
                           fill_value=np.nan).astype(np.float64)
            win_tf = src.window_transform(win)
            nodata = src.nodata

        nodata_mask = ~np.isfinite(arr)
        if nodata is not None:
            nodata_mask |= (arr == nodata)

        # Preserve NaN for cells that arrive as NoData (the resistance
        # raster uses NoData to mark IMPASSABLE pixels: highways,
        # large lakes, permanent fences). MCP_Geometric and
        # Circuitscape.jl both treat NaN as a wall, so LCPs and current
        # flow cannot cross. Finite cells are clipped to a numerically
        # sane range.
        valid = ~nodata_mask
        if valid.any():
            arr[valid] = np.clip(arr[valid], 1e-3, 1e6)
        arr[nodata_mask] = np.nan

        self.nodata_mask = nodata_mask
        return arr, win_tf, crs, arr.shape[0], arr.shape[1]

    def _apply_boar_resistance_penalty(self, arr):
        """Cubic penalty so wild boar treat high-R pixels as near-walls.

        Impassable cells (NaN) are PRESERVED as NaN throughout, so they
        remain walls to MCP_Geometric / Circuitscape downstream. The
        penalty is computed from the median / max of the FINITE cells
        only - otherwise NaN would contaminate the percentile statistics.
        """
        valid = (~self.nodata_mask if self.nodata_mask is not None
                 else np.isfinite(arr))
        finite_vals = arr[valid]
        if finite_vals.size == 0:
            return arr
        R_typ = float(np.median(finite_vals))
        if R_typ <= 0 or not np.isfinite(R_typ):
            R_typ = 1.0
        R_max = float(np.max(finite_vals))

        eff = arr * (arr / R_typ) ** 2
        # Clip ONLY finite cells; leave NaN walls intact.
        eff_valid = eff[valid]
        np.clip(eff_valid, 1e-3, 1e9, out=eff_valid)
        eff[valid] = eff_valid
        eff[~valid] = np.nan

        n_walls = int((~valid).sum())
        QgsMessageLog.logMessage(
            f"Boar penalty: R_typical={R_typ:.2f}, R_max={R_max:.2f}, "
            f"impassable cells (NaN walls): {n_walls}",
            LOG_TAG, Qgis.Info)
        return eff

    # -----------------------------------------------------------------
    # Origin disc rasterisation, with snap-to-passable fallback
    # -----------------------------------------------------------------
    def _origin_disc_mask(self, shape, win_tf, arr=None):
        """Flat-index array of the origin disc, restricted to passable cells.

        If `arr` is provided and the click landed on an impassable (NaN)
        cell or outside the raster, the origin is snapped to the nearest
        passable cell. `self.origin_xy` is updated to the snapped point
        and `self.origin_was_snapped` is set so the UI thread can flag it.
        """
        rows, cols = shape
        try:
            inv = ~win_tf
        except Exception:
            return np.array([], dtype=np.intp)
        c0, r0 = inv * self.origin_xy
        r0 = int(round(r0))
        c0 = int(round(c0))

        in_bounds = (0 <= r0 < rows) and (0 <= c0 < cols)
        passable_here = (in_bounds and arr is not None
                         and np.isfinite(arr[r0, c0]))

        # If the click is on a wall (or outside the data), snap.
        if arr is not None and not passable_here:
            snap = self._snap_to_passable(arr, r0, c0)
            if snap is None:
                return np.array([], dtype=np.intp)
            r0, c0 = snap
            x_snap, y_snap = rio_xy(win_tf, r0, c0)
            QgsMessageLog.logMessage(
                "Origin click was on an impassable cell. Snapped: "
                f"({self.origin_xy[0]:.1f}, {self.origin_xy[1]:.1f}) -> "
                f"({x_snap:.1f}, {y_snap:.1f}).",
                LOG_TAG, Qgis.Warning)
            self.origin_xy = (float(x_snap), float(y_snap))
            self.origin_was_snapped = True

        if not (0 <= r0 < rows and 0 <= c0 < cols):
            return np.array([], dtype=np.intp)

        radius = self.origin_radius_cells
        cells = []
        for dr in range(-radius, radius + 1):
            for dc in range(-radius, radius + 1):
                if dr * dr + dc * dc > radius * radius:
                    continue
                r, c = r0 + dr, c0 + dc
                if 0 <= r < rows and 0 <= c < cols:
                    # Skip any wall pixels that sneak into the disc.
                    if arr is None or np.isfinite(arr[r, c]):
                        cells.append(r * cols + c)
        return np.array(cells, dtype=np.intp)

    def _snap_to_passable(self, arr, r0, c0, max_radius_cells=100):
        """Find the (r, c) of the nearest passable cell to (r0, c0).

        Searches a square window of side 2*max_radius_cells (default
        ~2.5 km at 25 m/px) and returns the Euclidean-nearest cell with
        a finite resistance. Returns None if no such cell exists in the
        search window.
        """
        rows, cols = arr.shape
        # Clamp the search centre into the raster.
        rc = max(0, min(rows - 1, r0))
        cc = max(0, min(cols - 1, c0))
        if np.isfinite(arr[rc, cc]):
            return (rc, cc)

        r1 = max(0, rc - max_radius_cells)
        r2 = min(rows, rc + max_radius_cells + 1)
        c1 = max(0, cc - max_radius_cells)
        c2 = min(cols, cc + max_radius_cells + 1)
        sub = arr[r1:r2, c1:c2]
        passable = np.isfinite(sub)
        if not passable.any():
            return None
        rr, ccs = np.where(passable)
        # Distance from the (relative) click point.
        dr = rr - (rc - r1)
        dc = ccs - (cc - c1)
        d2 = dr * dr + dc * dc
        k = int(np.argmin(d2))
        return (int(rr[k] + r1), int(ccs[k] + c1))

    # -----------------------------------------------------------------
    def _single_source_dijkstra(self, arr, origin_flat_ids):
        """One Dijkstra sweep from the origin disc through the whole grid.

        IMPORTANT: skimage.graph.MCP_Geometric is implemented in C and
        causes an access-violation crash when the cost array contains
        NaN. It DOES handle np.inf correctly (treats those cells as
        impassable, just like NaN was supposed to behave). So we
        substitute NaN -> np.inf only for the Dijkstra call. The
        original arr (with NaN) is preserved for everything else.
        """
        from skimage.graph import MCP_Geometric
        rows, cols = arr.shape
        rr, cc = np.unravel_index(origin_flat_ids, (rows, cols))
        starts = list(zip(rr.tolist(), cc.tolist()))

        arr_for_mcp = np.where(np.isfinite(arr), arr, np.inf)
        mcp = MCP_Geometric(arr_for_mcp, fully_connected=True)
        cost_grid, _ = mcp.find_costs(starts)
        return cost_grid, mcp

    # -----------------------------------------------------------------
    def _auto_aoi_mask(self, cost_grid, origin_flat_ids, shape):
        """Cells reachable from origin with cost below an auto threshold.

        Threshold = 4 x median(finite costs). With the cubic boar penalty
        applied, this naturally cuts off areas behind barriers without a
        user-defined radius.
        """
        rows, cols = shape
        finite = cost_grid[np.isfinite(cost_grid)]
        if finite.size == 0:
            aoi = np.zeros(shape, dtype=bool)
        else:
            median = float(np.median(finite))
            threshold = max(median * 4.0,
                            float(np.percentile(finite, 60)) * 1.5)
            aoi = np.isfinite(cost_grid) & (cost_grid <= threshold)
        # Always include the origin disc.
        rr, cc = np.unravel_index(origin_flat_ids, (rows, cols))
        aoi[rr, cc] = True
        return aoi

    # -----------------------------------------------------------------
    # AOI boundary cells - used as the Circuit-theory sink.
    # -----------------------------------------------------------------
    def _aoi_boundary_mask(self):
        aoi = self.aoi_mask
        if aoi is None:
            return None
        not_aoi = ~aoi
        up    = np.zeros_like(aoi); up[:-1]    = not_aoi[1:]
        down  = np.zeros_like(aoi); down[1:]   = not_aoi[:-1]
        left  = np.zeros_like(aoi); left[:, :-1]  = not_aoi[:, 1:]
        right = np.zeros_like(aoi); right[:, 1:]  = not_aoi[:, :-1]
        return aoi & (up | down | left | right)

    # -----------------------------------------------------------------
    # LCPs to FOREST destinations (auto-detected from the resistance raster).
    #
    # Method:
    #   1. Threshold the in-AOI valid cells at the 20th percentile to
    #      isolate "forest" pixels.
    #   2. Find connected components (scipy.ndimage.label, 8-connected).
    #   3. Skip the patch(es) overlapping the origin disc.
    #   4. Pick the LARGEST N remaining patches; anchor each by its
    #      centroid (or nearest in-patch cell if the centroid lies
    #      outside a concave patch).
    #   5. Traceback the MCP from each anchor.
    #   6. Edge traffic = number of LCPs sharing the edge (+1 per LCP).
    #   7. Merge consecutive edges with equal traffic into single
    #      line features.
    # -----------------------------------------------------------------
    def _build_lcps_to_forests(self, arr, win_tf, mcp,
                               n_max_targets=12,
                               forest_quantile=0.20,
                               min_patch_cells=5):
        try:
            from scipy.ndimage import label
        except ImportError:
            QgsMessageLog.logMessage(
                "scipy.ndimage missing; LCPs skipped.",
                LOG_TAG, Qgis.Warning)
            return

        rows, cols = arr.shape

        # --- Build the forest mask -----------------------------------
        in_play = np.ones(arr.shape, dtype=bool)
        if self.aoi_mask is not None:
            in_play &= self.aoi_mask
        if self.nodata_mask is not None:
            in_play &= ~self.nodata_mask
        valid_values = arr[in_play]
        if valid_values.size == 0:
            QgsMessageLog.logMessage(
                "No valid in-AOI cells; LCPs skipped.",
                LOG_TAG, Qgis.Warning)
            return
        threshold = float(np.percentile(valid_values,
                                        100.0 * forest_quantile))
        forest_mask = (arr <= threshold) & in_play

        if not forest_mask.any():
            QgsMessageLog.logMessage(
                "No forest pixels under threshold; LCPs skipped.",
                LOG_TAG, Qgis.Warning)
            return

        # --- Connected components (8-connected, like our Dijkstra) ---
        labeled, n_features = label(forest_mask,
                                    structure=np.ones((3, 3), dtype=int))
        if n_features == 0:
            return

        # --- Identify the patch(es) overlapping the origin disc ------
        origin_ids = self._origin_disc_mask((rows, cols), win_tf, arr=arr)
        origin_patch_ids = set()
        if origin_ids.size:
            orr, occ = np.unravel_index(origin_ids, (rows, cols))
            lab_at_origin = labeled[orr, occ]
            origin_patch_ids = {int(v) for v in lab_at_origin if v > 0}

        # --- Anchor cell per (non-origin) patch ----------------------
        patches = []   # list of (anchor_row, anchor_col, size)
        for cid in range(1, n_features + 1):
            if cid in origin_patch_ids:
                continue
            rr_p, cc_p = np.where(labeled == cid)
            size = rr_p.size
            if size < min_patch_cells:
                continue
            cr = int(round(float(rr_p.mean())))
            cc = int(round(float(cc_p.mean())))
            if labeled[cr, cc] != cid:
                # Concave patch: snap centroid to the closest in-patch cell.
                d2 = (rr_p - cr) ** 2 + (cc_p - cc) ** 2
                k = int(np.argmin(d2))
                cr, cc = int(rr_p[k]), int(cc_p[k])
            patches.append((cr, cc, size))

        if not patches:
            QgsMessageLog.logMessage(
                "Only the origin's own forest is in the AOI; no LCPs.",
                LOG_TAG, Qgis.Info)
            return

        # Largest forests first; cap at n_max_targets.
        patches.sort(key=lambda p: -p[2])
        patches = patches[:int(n_max_targets)]
        QgsMessageLog.logMessage(
            f"LCP forest destinations: {len(patches)} patches "
            f"(sizes: {[p[2] for p in patches[:5]]}...).",
            LOG_TAG, Qgis.Info)

        # Remember anchors for the visualisation layer (in world coords).
        self.forest_anchors = [
            (rio_xy(win_tf, r, c), size) for (r, c, size) in patches
        ]

        # --- Traceback each anchor -----------------------------------
        paths = []
        for r, c, _size in patches:
            try:
                indices = mcp.traceback((r, c))
            except Exception:
                continue
            cells = [(int(ri), int(ci)) for ri, ci in indices]
            if len(cells) >= 2:
                paths.append(cells)
        self.n_lcps = len(paths)
        if not paths:
            return

        # --- Merge into uniform-traffic polylines --------------------
        # Each shared edge accumulates +1 per LCP that uses it.
        from collections import defaultdict
        edge_count = defaultdict(int)
        for cells in paths:
            for a, b in zip(cells, cells[1:]):
                edge = (a, b) if a < b else (b, a)
                edge_count[edge] += 1

        visited = set()
        runs = []
        for cells in paths:
            seg_cells = []
            seg_count = None
            for a, b in zip(cells, cells[1:]):
                edge = (a, b) if a < b else (b, a)
                if edge in visited:
                    if seg_count is not None and len(seg_cells) >= 2:
                        runs.append((seg_cells, seg_count))
                    seg_cells = []
                    seg_count = None
                    continue
                visited.add(edge)
                c = edge_count[edge]
                if seg_count is None:
                    seg_cells = [a, b]
                    seg_count = c
                elif c == seg_count:
                    seg_cells.append(b)
                else:
                    runs.append((seg_cells, seg_count))
                    seg_cells = [a, b]
                    seg_count = c
            if seg_count is not None and len(seg_cells) >= 2:
                runs.append((seg_cells, seg_count))

        for cells, count in runs:
            pts = [QgsPointXY(*rio_xy(win_tf, r, c)) for r, c in cells]
            self.merged_lcp_polylines.append({
                "points":  pts,
                "traffic": int(count),
            })

    # =================================================================
    # Output rasters
    # =================================================================
    def _build_cost_raster(self, cost_grid, win_tf):
        """Optional: raw cost-from-origin (Dijkstra cost grid)."""
        return self._mask_crop_write_raster(
            cost_grid.astype(np.float32), win_tf, prefix="cost")

    def _build_risk_raster(self, cost_grid, win_tf):
        """Continuous infection-risk raster: exp(-cost / median_cost)."""
        rows, cols = cost_grid.shape
        in_aoi = cost_grid[self.aoi_mask] if self.aoi_mask is not None else cost_grid
        finite = in_aoi[np.isfinite(in_aoi)]
        scale = float(np.median(finite)) if finite.size else 1.0
        if scale <= 0:
            scale = 1.0
        risk = np.exp(-cost_grid / scale)
        risk[~np.isfinite(cost_grid)] = np.nan
        return self._mask_crop_write_raster(
            risk.astype(np.float32), win_tf, prefix="risk")

    # -----------------------------------------------------------------
    # Pinchpoint raster: Circuit theory with AOI boundary as sink.
    # -----------------------------------------------------------------
    def _try_circuitscape_jl(self, arr, win_tf, origin_flat_ids):
        rows, cols = arr.shape
        source_mask = np.zeros((rows, cols), dtype=bool)
        srr, scc = np.unravel_index(origin_flat_ids, (rows, cols))
        source_mask[srr, scc] = True

        sink_mask = self._aoi_boundary_mask()
        if sink_mask is None or not sink_mask.any():
            return False

        cs_resistance = arr.astype(np.float64, copy=True)
        if self.nodata_mask is not None:
            cs_resistance[self.nodata_mask] = np.nan
        if self.aoi_mask is not None:
            cs_resistance[~self.aoi_mask] = np.nan

        def _log(msg):
            QgsMessageLog.logMessage(msg, LOG_TAG, Qgis.Info)

        current = run_circuitscape_jl(
            cs_resistance, win_tf,
            source_mask=source_mask,
            sink_mask=sink_mask,
            log=_log,
        )

        finite = current[np.isfinite(current) & (current > 0)]
        floor = max(float(np.percentile(finite, 1)), 1e-12) \
                if finite.size else 1e-12
        with np.errstate(divide="ignore", invalid="ignore"):
            log_curr = np.log10(np.maximum(current, floor))
        log_curr[~np.isfinite(current)] = np.nan
        self.current_raster_path = self._mask_crop_write_raster(
            log_curr.astype(np.float32), win_tf, prefix="pinchpoints_cs")
        return True

    def _single_source_current_scipy(self, arr, win_tf, crs, origin_flat_ids):
        """scipy fallback: assemble L, solve once for origin disc vs AOI border.

        NaN (impassable) cells would create zero-degree isolated nodes
        and make the Laplacian singular -> spsolve crashes. We
        substitute NaN with 1e10 here so the system stays non-singular;
        conductance through those cells becomes ~1e-10, which is
        effectively a wall (same physical behaviour, no numerical
        explosion). Circuitscape.jl handles NaN natively and is the
        preferred path when Julia is available.
        """
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import spsolve

        rows, cols = arr.shape
        n = rows * cols
        R_flat = arr.flatten()
        # Replace NaN walls with finite-but-huge resistance.
        R_flat = np.where(np.isfinite(R_flat), R_flat, 1e10)

        idx = np.arange(n).reshape(rows, cols)
        h_i = idx[:, :-1].flatten();  h_j = idx[:, 1:].flatten()
        v_i = idx[:-1, :].flatten();  v_j = idx[1:, :].flatten()
        e_i = np.concatenate([h_i, v_i])
        e_j = np.concatenate([h_j, v_j])

        R_edge = 0.5 * (R_flat[e_i] + R_flat[e_j])
        # All edges are now finite by construction; no need to filter.
        g_edge = 1.0 / np.maximum(R_edge, 1e-9)

        data = np.concatenate([-g_edge, -g_edge,  g_edge,  g_edge])
        row  = np.concatenate([   e_i,    e_j,    e_i,    e_j])
        col  = np.concatenate([   e_j,    e_i,    e_i,    e_j])
        L = coo_matrix((data, (row, col)), shape=(n, n)).tocsc()

        # Source: origin disc.
        bvec = np.zeros(n)
        bvec[origin_flat_ids] += 1.0 / origin_flat_ids.size

        # Sink: AOI boundary cells.
        sink_mask = self._aoi_boundary_mask()
        if sink_mask is None or not sink_mask.any():
            return
        sink_ids = np.flatnonzero(sink_mask.flatten())
        bvec[sink_ids] -= 1.0 / sink_ids.size

        ground = int(sink_ids[0])
        keep = np.ones(n, dtype=bool)
        keep[ground] = False
        L_red = L[keep][:, keep]
        b_red = bvec[keep]

        v_red = spsolve(L_red, b_red)
        v = np.zeros(n)
        v[keep] = v_red

        edge_I = g_edge * (v[e_i] - v[e_j])
        node_current = np.zeros(n)
        np.add.at(node_current, e_i, np.abs(edge_I))
        np.add.at(node_current, e_j, np.abs(edge_I))
        node_current *= 0.5
        node_current = node_current.reshape(rows, cols).astype(np.float64)

        finite = node_current[node_current > 0]
        floor = max(float(np.percentile(finite, 1)), 1e-12) \
                if finite.size else 1e-12
        with np.errstate(divide="ignore", invalid="ignore"):
            log_curr = np.log10(np.maximum(node_current, floor))
        log_curr = log_curr.astype(np.float32)
        self.current_raster_path = self._mask_crop_write_raster(
            log_curr, win_tf, prefix="pinchpoints")

    # -----------------------------------------------------------------
    # Random walks from the origin disc.
    #
    # Biological / numerical correctness:
    #   - 8-connected, but diagonal moves cost sqrt(2) more than
    #     orthogonal ones. We weight by 1 / (R_neighbour * step_distance)
    #     so the walker doesn't artificially prefer diagonals.
    #   - NoData neighbours are excluded outright (a boar never steps
    #     into a cell with no resistance data).
    #   - Walk terminates when leaving the AOI - that's the realistic
    #     boundary of boar dispersal.
    # The output raster sits on the same pixel grid as the resistance
    # raster (window_transform preserves alignment), so overlays in QGIS
    # line up cell-for-cell.
    # -----------------------------------------------------------------
    def _random_walks(self, arr, win_tf, crs, origin_flat_ids,
                      n_walks=200, max_steps=2000):
        rows, cols = arr.shape
        rr, cc = np.unravel_index(origin_flat_ids, (rows, cols))
        starts = np.column_stack([rr, cc])

        visit = np.zeros((rows, cols), dtype=np.int32)
        rng = np.random.default_rng()

        offsets = np.array([
            (-1, -1), (-1, 0), (-1, 1),
            ( 0, -1),          ( 0, 1),
            ( 1, -1), ( 1, 0), ( 1, 1),
        ], dtype=np.intp)
        # Step distance per offset: orthogonal = 1, diagonal = sqrt(2).
        SQ2 = float(np.sqrt(2.0))
        step_dist = np.array([
            SQ2, 1.0, SQ2,
            1.0,      1.0,
            SQ2, 1.0, SQ2,
        ], dtype=np.float64)

        aoi    = self.aoi_mask
        nodata = self.nodata_mask

        for _ in range(n_walks):
            if self.isCanceled():
                return
            r, c = starts[rng.integers(0, len(starts))]
            r, c = int(r), int(c)
            for _ in range(max_steps):
                visit[r, c] += 1
                # Stop once outside the AOI.
                if aoi is not None and not aoi[r, c]:
                    break

                nrs = r + offsets[:, 0]
                ncs = c + offsets[:, 1]
                ok = (nrs >= 0) & (nrs < rows) \
                   & (ncs >= 0) & (ncs < cols)

                # Exclude NoData neighbours (boars don't walk into "no data").
                if nodata is not None and ok.any():
                    safe_r = np.clip(nrs, 0, rows - 1)
                    safe_c = np.clip(ncs, 0, cols - 1)
                    ok = ok & ~nodata[safe_r, safe_c]

                if not ok.any():
                    break

                nrs_ok = nrs[ok]
                ncs_ok = ncs[ok]
                dists  = step_dist[ok]
                # Step cost = R_neighbour * step_distance; pick neighbour
                # with probability inversely proportional to that cost.
                w = 1.0 / (arr[nrs_ok, ncs_ok] * dists)
                w /= w.sum()
                k = rng.choice(len(nrs_ok), p=w)
                r, c = int(nrs_ok[k]), int(ncs_ok[k])

        if visit.any():
            grid = visit.astype(np.float32)
            grid[grid == 0] = np.nan
            self.walk_raster_path = self._mask_crop_write_raster(
                grid, win_tf, prefix="random_walk")

    # =================================================================
    # Mask + crop + write helper
    # =================================================================
    def _mask_crop_write_raster(self, grid, win_tf, prefix):
        """Apply AOI + NoData mask, crop tightly, save GeoTIFF."""
        grid = grid.astype(np.float32, copy=True)

        if self.aoi_mask is not None and self.aoi_mask.shape == grid.shape:
            inside = self.aoi_mask
        else:
            inside = np.ones(grid.shape, dtype=bool)
        if self.nodata_mask is not None \
                and self.nodata_mask.shape == grid.shape:
            inside = inside & ~self.nodata_mask

        grid[~inside] = np.nan

        v_rows = np.where(inside.any(axis=1))[0]
        v_cols = np.where(inside.any(axis=0))[0]
        if not v_rows.size or not v_cols.size:
            return None
        r0, r1 = int(v_rows[0]), int(v_rows[-1]) + 1
        c0, c1 = int(v_cols[0]), int(v_cols[-1]) + 1
        grid = grid[r0:r1, c0:c1]
        out_tf = window_transform(Window(c0, r0, c1 - c0, r1 - r0), win_tf)
        out_h, out_w = grid.shape

        out_arr = np.where(np.isnan(grid), -9999.0, grid).astype(np.float32)
        out_path = os.path.join(
            tempfile.gettempdir(),
            f"wildboar_{prefix}_{os.getpid()}.tif")
        from rasterio.crs import CRS
        try:
            crs = CRS.from_wkt(self.crs_wkt)
        except Exception:
            crs = None
        with rasterio.open(
            out_path, "w",
            driver="GTiff",
            height=out_h, width=out_w, count=1,
            dtype="float32",
            crs=crs, transform=out_tf,
            nodata=-9999,
        ) as dst:
            dst.write(out_arr, 1)
        return out_path

    # =================================================================
    # UI thread: finalisation
    # =================================================================
    def finished(self, ok):
        if not ok:
            msg = (f"Failed: {self.exception}" if self.exception
                   else "Cancelled.")
            QMessageBox.critical(None, "Wildboar Connectivity", msg)
            return

        # If the origin was snapped from an impassable click, show the
        # actual analysis origin so the user can see what happened.
        if self.origin_was_snapped:
            self._add_snapped_origin_layer()

        # LCP corridors as merged lines (one feature per uniform-traffic run)
        if self.merged_lcp_polylines:
            self._add_lcp_lines_layer()
        # Visual anchors showing what the LCPs were aimed at.
        if self.forest_anchors:
            self._add_forest_anchors_layer()

        # Pinchpoint raster (the corridor map)
        if self.current_raster_path:
            r = QgsRasterLayer(self.current_raster_path,
                               "Dispersal corridor probability (Circuit theory, log10)")
            if r.isValid():
                self._apply_singleband_ramp(r, "Magma")
                add_wildboar_layer(r, ZOrder.PINCHPOINTS)

        # Continuous infection-risk raster
        if self.risk_raster_path:
            r = QgsRasterLayer(self.risk_raster_path,
                               "Cumulative dispersal probability (resistant kernel)")
            if r.isValid():
                self._apply_risk_ramp(r)
                add_wildboar_layer(r, ZOrder.SELECTED_HABITATS + 1)

        QgsMessageLog.logMessage(
            "Done. lcps={l}({s} segments) pinch={p} risk={r}".format(
                l=self.n_lcps,
                s=len(self.merged_lcp_polylines),
                p="yes" if self.current_raster_path else "no",
                r="yes" if self.risk_raster_path else "no"),
            LOG_TAG, Qgis.Success)

    # -----------------------------------------------------------------
    def _add_lcp_lines_layer(self):
        """Merged LCP polylines styled by traffic count in a BLUE palette.

        The blue family was chosen deliberately to remain visible against
        both the green-yellow-red infection-risk classes and the
        black-purple-pink-yellow Magma pinchpoint raster. Neither
        background contains pure blue, so the corridors never "blend in".
        Three classes:
            1 path        -> light sky blue, thin
            mid traffic   -> vivid blue,     medium
            top traffic   -> deep navy,      thick (corridor backbone)
        """
        from qgis.core import (
            QgsGraduatedSymbolRenderer,
            QgsRendererRange,
            QgsSymbol,
        )
        polys = self.merged_lcp_polylines
        if not polys:
            return
        layer = QgsVectorLayer(
            f"LineString?crs={self.crs_wkt}",
            f"Example dispersal trajectories ({len(polys)} segments, "
            f"{self.n_lcps} target forests) — illustrative",
            "memory")
        prov = layer.dataProvider()
        prov.addAttributes([QgsField("traffic", QVariant.Int)])
        layer.updateFields()

        for p in polys:
            f = QgsFeature()
            f.setGeometry(QgsGeometry.fromPolylineXY(p["points"]))
            f.setAttributes([int(p["traffic"])])
            prov.addFeature(f)
        layer.updateExtents()

        traffics = sorted(int(p["traffic"]) for p in polys)
        max_t = traffics[-1] if traffics else 1
        if max_t <= 1:
            low_hi, med_hi = 1.5, 1.5
        else:
            q33 = traffics[len(traffics) // 3]
            q66 = traffics[(2 * len(traffics)) // 3]
            low_hi = max(1.5, q33 + 0.5)
            med_hi = max(low_hi + 0.5, q66 + 0.5)

        def line_sym(color, width):
            s = QgsSymbol.defaultSymbol(layer.geometryType())
            s.setColor(QColor(color))
            s.setWidth(width)
            return s

        upper = max_t + 0.5
        ranges = [
            QgsRendererRange(0.5, low_hi,
                             line_sym("#85C1E9", 0.6),
                             "1 path (detour)"),
            QgsRendererRange(low_hi, med_hi,
                             line_sym("#2874A6", 1.6),
                             f"{int(low_hi+0.5)}-{int(med_hi-0.5)} paths"),
            QgsRendererRange(med_hi, upper,
                             line_sym("#154360", 3.2),
                             f">= {int(med_hi+0.5)} paths (corridor!)"),
        ]
        renderer = QgsGraduatedSymbolRenderer("traffic", ranges)
        layer.setRenderer(renderer)
        add_wildboar_layer(layer, ZOrder.LCP_TRAFFIC)

    # -----------------------------------------------------------------
    def _add_snapped_origin_layer(self):
        """Marker showing where the analysis origin actually sits after
        the impassable-click snap. The original click point is still
        displayed by main_plugin as the red 'Outbreak origin' marker;
        this orange marker shows the snapped position so the user can
        see both."""
        layer = QgsVectorLayer(
            f"Point?crs={self.crs_wkt}",
            "Outbreak origin (snapped to passable cell)",
            "memory")
        prov = layer.dataProvider()
        f = QgsFeature()
        f.setGeometry(QgsGeometry.fromPointXY(
            QgsPointXY(self.origin_xy[0], self.origin_xy[1])))
        prov.addFeature(f)
        layer.updateExtents()
        sym = layer.renderer().symbol()
        sym.setColor(QColor("#f39c12"))   # orange - distinct from the red click
        sym.setSize(5.0)
        add_wildboar_layer(layer, ZOrder.CENTRE)

    # -----------------------------------------------------------------
    def _add_forest_anchors_layer(self):
        """Show the auto-detected forest destinations as labelled dots."""
        layer = QgsVectorLayer(
            f"Point?crs={self.crs_wkt}",
            f"Forest destinations ({len(self.forest_anchors)})",
            "memory")
        prov = layer.dataProvider()
        prov.addAttributes([QgsField("size_cells", QVariant.Int)])
        layer.updateFields()
        for (x, y), size in self.forest_anchors:
            f = QgsFeature()
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
            f.setAttributes([int(size)])
            prov.addFeature(f)
        layer.updateExtents()
        sym = layer.renderer().symbol()
        sym.setColor(QColor("#27ae60"))     # forest green
        sym.setSize(3.5)
        add_wildboar_layer(layer, ZOrder.SELECTED_HABITATS)

    # =================================================================
    # Renderers
    # =================================================================
    @staticmethod
    def _apply_singleband_ramp(rlayer, ramp_name="Magma"):
        from qgis.core import (
            QgsSingleBandPseudoColorRenderer,
            QgsColorRampShader,
            QgsRasterShader,
            QgsStyle,
        )
        try:
            with rasterio.open(rlayer.source()) as ds:
                data = ds.read(1, masked=True)
            valid = data.compressed()
            if valid.size:
                vmin = float(np.nanpercentile(valid, 2))
                vmax = float(np.nanpercentile(valid, 98))
            else:
                raise ValueError("no valid pixels")
        except Exception:
            stats = rlayer.dataProvider().bandStatistics(1)
            vmin, vmax = stats.minimumValue, stats.maximumValue
        if not (np.isfinite(vmin) and np.isfinite(vmax)) or vmax <= vmin:
            vmin, vmax = 0.0, 1.0

        style = QgsStyle.defaultStyle()
        ramp = style.colorRamp(ramp_name) or style.colorRamp("Viridis")
        shader_fn = QgsColorRampShader(vmin, vmax, ramp,
                                       QgsColorRampShader.Interpolated)
        shader_fn.classifyColorRamp()
        shader = QgsRasterShader()
        shader.setRasterShaderFunction(shader_fn)
        renderer = QgsSingleBandPseudoColorRenderer(
            rlayer.dataProvider(), 1, shader)
        rlayer.setRenderer(renderer)
        rlayer.triggerRepaint()

    @staticmethod
    def _apply_risk_ramp(rlayer):
        """Discrete three-band infection-risk classification.

        Cutoffs are derived from the boar-dispersal interpretation of
        the risk function   risk = exp(-cost / median_cost) :

            HIGH    risk > 0.50   <=>   cost < 0.69 x median
                    Boars routinely reach this cell during normal
                    movement. Almost guaranteed ASF exposure.

            MEDIUM  0.15 < risk <= 0.50  <=>  0.69 x median < cost < 1.9 x median
                    Reachable with extended movement (dispersal-scale).
                    Plausible ASF exposure within months.

            LOW     risk <= 0.15  <=>   cost > 1.9 x median
                    Within the AOI but at its far edge - unlikely to be
                    reached by a single boar.

        The cutoffs correspond to e^-0.69 and e^-1.9 - half-life and
        decay-to-15%, both standard breakpoints for exponential decay
        risk functions in spatial epidemiology.
        """
        from qgis.core import (
            QgsColorRampShader,
            QgsRasterShader,
            QgsSingleBandPseudoColorRenderer,
        )
        shader_fn = QgsColorRampShader()
        shader_fn.setColorRampType(QgsColorRampShader.Discrete)
        shader_fn.setColorRampItemList([
            QgsColorRampShader.ColorRampItem(
                0.15, QColor("#2ecc71"), "Low risk (risk ≤ 0.15)"),
            QgsColorRampShader.ColorRampItem(
                0.50, QColor("#f1c40f"), "Medium risk (0.15 < risk ≤ 0.50)"),
            QgsColorRampShader.ColorRampItem(
                1.01, QColor("#e74c3c"), "High risk (risk > 0.50)"),
        ])
        shader = QgsRasterShader()
        shader.setRasterShaderFunction(shader_fn)
        renderer = QgsSingleBandPseudoColorRenderer(
            rlayer.dataProvider(), 1, shader)
        renderer.setOpacity(0.65)
        rlayer.setRenderer(renderer)
        rlayer.triggerRepaint()
