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
            origin_cells = self._origin_disc_mask(arr.shape, win_tf)
            if origin_cells.size == 0:
                raise RuntimeError(
                    "Outbreak point falls outside the resistance raster.")
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

            # ---- Cost-from-origin raster (optional) ------------------
            if self.options.get("cost", False):
                self.cost_raster_path = self._build_cost_raster(
                    cost_grid, win_tf)
            self.setProgress(33)
            if self.isCanceled():
                return False

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

            # ---- Optional: random walks ------------------------------
            if self.options.get("random_walk", False):
                self._random_walks(
                    arr, win_tf, crs, origin_cells,
                    n_walks=int(self.options.get("n_walks", 200)))

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
        arr[nodata_mask] = 1e6
        np.clip(arr, 1e-3, 1e6, out=arr)

        self.nodata_mask = nodata_mask
        return arr, win_tf, crs, arr.shape[0], arr.shape[1]

    def _apply_boar_resistance_penalty(self, arr):
        """Cubic penalty so wild boar treat high-R pixels as near-walls."""
        if self.nodata_mask is not None and self.nodata_mask.any():
            valid = arr[~self.nodata_mask]
        else:
            valid = arr
        if valid.size == 0:
            return arr
        R_typ = float(np.median(valid))
        if R_typ <= 0 or not np.isfinite(R_typ):
            R_typ = 1.0
        R_max = float(np.max(valid))

        eff = arr * (arr / R_typ) ** 2
        np.clip(eff, 1e-3, 1e9, out=eff)
        if self.nodata_mask is not None:
            eff[self.nodata_mask] = 1e9

        QgsMessageLog.logMessage(
            f"Boar penalty: R_typical={R_typ:.2f}, R_max={R_max:.2f}, "
            f"max-cell costs {(R_max ** 2 / R_typ ** 2):.0f}x a typical cell.",
            LOG_TAG, Qgis.Info)
        return eff

    # -----------------------------------------------------------------
    # Origin disc rasterisation
    # -----------------------------------------------------------------
    def _origin_disc_mask(self, shape, win_tf):
        """Return array of flat indices for the origin disc on the window."""
        rows, cols = shape
        try:
            inv = ~win_tf
        except Exception:
            return np.array([], dtype=np.intp)
        c0, r0 = inv * self.origin_xy
        r0 = int(round(r0))
        c0 = int(round(c0))
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
                    cells.append(r * cols + c)
        return np.array(cells, dtype=np.intp)

    # -----------------------------------------------------------------
    def _single_source_dijkstra(self, arr, origin_flat_ids):
        from skimage.graph import MCP_Geometric
        rows, cols = arr.shape
        rr, cc = np.unravel_index(origin_flat_ids, (rows, cols))
        starts = list(zip(rr.tolist(), cc.tolist()))
        mcp = MCP_Geometric(arr, fully_connected=True)
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
    # Without core-habitat polygons, the destinations must come from the
    # raster itself. Boars disperse toward attractive habitat, which on a
    # boar-resistance surface is the LOW-resistance cells (forest). We:
    #
    #   1. Threshold the in-AOI valid cells at the 20th percentile to
    #      isolate "forest" pixels.
    #   2. Find connected components (scipy.ndimage.label, 8-connected).
    #   3. Skip the patch(es) overlapping the origin disc - those are
    #      the boar's own forest and "leading to itself" makes no sense.
    #   4. Pick the largest N remaining patches and snap an anchor cell
    #      to each (centroid, or nearest in-patch cell if the centroid
    #      lies outside a concave patch).
    #   5. Traceback through the MCP from each anchor.
    #   6. Merge overlapping LCP segments into uniform-traffic runs.
    #
    # Output also includes a small "Forest destinations" point layer so
    # the user can see WHERE the LCPs were aimed.
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
        origin_ids = self._origin_disc_mask((rows, cols), win_tf)
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

        # --- Traceback each anchor --------------------------------------
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

        # --- Merge into uniform-traffic polylines ---------------------
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
        """scipy fallback: assemble L, solve once for origin disc vs AOI border."""
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import spsolve

        rows, cols = arr.shape
        n = rows * cols
        R_flat = arr.flatten()

        idx = np.arange(n).reshape(rows, cols)
        h_i = idx[:, :-1].flatten();  h_j = idx[:, 1:].flatten()
        v_i = idx[:-1, :].flatten();  v_j = idx[1:, :].flatten()
        e_i = np.concatenate([h_i, v_i])
        e_j = np.concatenate([h_j, v_j])

        R_edge = 0.5 * (R_flat[e_i] + R_flat[e_j])
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

        # LCP corridors as merged lines (one feature per uniform-traffic run)
        if self.merged_lcp_polylines:
            self._add_lcp_lines_layer()
        # Visual anchors showing what the LCPs were aimed at.
        if self.forest_anchors:
            self._add_forest_anchors_layer()

        # Pinchpoint raster (the corridor map)
        if self.current_raster_path:
            r = QgsRasterLayer(self.current_raster_path,
                               "Pinchpoints / corridors (log10 current)")
            if r.isValid():
                self._apply_singleband_ramp(r, "Magma")
                add_wildboar_layer(r, ZOrder.PINCHPOINTS)

        # Continuous infection-risk raster
        if self.risk_raster_path:
            r = QgsRasterLayer(self.risk_raster_path,
                               "Infection risk (continuous)")
            if r.isValid():
                self._apply_risk_ramp(r)
                add_wildboar_layer(r, ZOrder.SELECTED_HABITATS + 1)

        # Optional cost-from-origin raster
        if self.cost_raster_path:
            r = QgsRasterLayer(self.cost_raster_path,
                               "Cost-from-origin (boar distance surface)")
            if r.isValid():
                self._apply_singleband_ramp(r, "Viridis")
                add_wildboar_layer(r, ZOrder.LCP_TRAFFIC)

        # Optional random walk density
        if self.walk_raster_path:
            r = QgsRasterLayer(
                self.walk_raster_path,
                f"Random walk density (n={int(self.options.get('n_walks', 200))})")
            if r.isValid():
                self._apply_singleband_ramp(r, "Viridis")
                add_wildboar_layer(r, ZOrder.GRAPH_EDGES)

        QgsMessageLog.logMessage(
            "Done. lcps={l}({s} segments) pinch={p} risk={r} "
            "cost={c} walks={w}".format(
                l=self.n_lcps,
                s=len(self.merged_lcp_polylines),
                p="yes" if self.current_raster_path else "no",
                r="yes" if self.risk_raster_path else "no",
                c="yes" if self.cost_raster_path else "no",
                w="yes" if self.walk_raster_path else "no"),
            LOG_TAG, Qgis.Success)

    # -----------------------------------------------------------------
    def _add_lcp_lines_layer(self):
        """Render merged LCP polylines, styled by traffic count.

        Three classes (thin grey -> orange -> thick dark red) make the
        corridor backbone visually obvious without further user styling.
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
            f"LCP corridors (merged, {len(polys)} segments, "
            f"{self.n_lcps} paths)",
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
                             line_sym("#95a5a6", 0.4),
                             "1 path (detour)"),
            QgsRendererRange(low_hi, med_hi,
                             line_sym("#e67e22", 1.2),
                             f"{int(low_hi+0.5)}-{int(med_hi-0.5)} paths"),
            QgsRendererRange(med_hi, upper,
                             line_sym("#c0392b", 2.6),
                             f">= {int(med_hi+0.5)} paths (corridor!)"),
        ]
        renderer = QgsGraduatedSymbolRenderer("traffic", ranges)
        layer.setRenderer(renderer)
        add_wildboar_layer(layer, ZOrder.LCP_TRAFFIC)

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
        from qgis.core import (
            QgsColorRampShader,
            QgsGradientColorRamp,
            QgsGradientStop,
            QgsRasterShader,
            QgsSingleBandPseudoColorRenderer,
        )
        low_color  = QColor("#2ecc71")
        mid_color  = QColor("#f1c40f")
        high_color = QColor("#e74c3c")
        ramp = QgsGradientColorRamp(low_color, high_color)
        ramp.setStops([QgsGradientStop(0.5, mid_color)])
        shader_fn = QgsColorRampShader(
            0.0, 1.0, ramp, QgsColorRampShader.Interpolated)
        shader_fn.classifyColorRamp()
        shader = QgsRasterShader()
        shader.setRasterShaderFunction(shader_fn)
        renderer = QgsSingleBandPseudoColorRenderer(
            rlayer.dataProvider(), 1, shader)
        renderer.setOpacity(0.75)
        rlayer.setRenderer(renderer)
        rlayer.triggerRepaint()
