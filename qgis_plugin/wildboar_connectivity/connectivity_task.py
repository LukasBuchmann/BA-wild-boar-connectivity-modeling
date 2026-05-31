# -*- coding: utf-8 -*-
"""
Background QgsTask for the ASF Wildboar Connectivity plugin.

WORKFLOW:

    Outbreak point  ->  small source disc (_point_source_cells)
                    ->  Dijkstra from source disc
                    ->  auto AOI  (all cells with finite cost)
                    ->  infection-risk raster     exp(-cost/D_cost)
                    ->  LCP corridors             to nearest low-resistance clusters
                    ->  pinchpoint raster         Circuitscape.jl or scipy fallback
                    ->  BCRW random-walk density  (optional)
                    ->  iSSF-IBMM contamination   (optional)
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


# ---------------------------------------------------------------------
# Biological calibration of the infection-risk decay kernel.
#
# risk(d) = exp(-d / D)  where d is cost-equivalent distance through
# typical-resistance terrain.
#
# D = 4 km, the consensus mean natal dispersal distance for European
# wild boar across the published telemetry / mark-recapture record:
#
#   Truve and Lemel (2003)  Wildlife Biology  9(4): 51 to 57.
#   Keuling et al. (2010)   Eur J Wildl Res  56(2): 159 to 167.
#   Prevot and Licoppe (2013) Eur J Wildl Res 59(6): 795 to 803.
#   Morelle, Lehaire, Lejeune (2015) Mammal Review 45: 15 to 29.
#
# The resulting risk bands map to the EFSA operational zones (2018
# Scientific Opinion on ASF in wild boar):
#   risk > 0.50 (within ~2.8 km)  -> EFSA protection zone (3 km).
#   risk 0.15 to 0.50 (~2.8 to ~7.6 km) -> EFSA surveillance zone (10 km).
#   risk <= 0.15  (> ~7.6 km)     -> outside the surveillance zone.
# ---------------------------------------------------------------------
ASF_DISPERSAL_KM = 4.0

# =================================================================
# iSSF-IBMM parameters
# Movement kernel fitted to WildBoar_InMatrix_Steps.csv 
# Updated with final iSSF coefficients (log_sl = 0.153, cos_ta = 0.564)
# =================================================================
IBMM_SL_GAMMA_SHAPE    = 2.65258  # tentative (2.5) + beta_log_sl
IBMM_SL_GAMMA_SCALE_M  = 160.0    # step-length scale (metres)
IBMM_KAPPA             = 1.06417  # tentative (0.5) + beta_cos_ta 
IBMM_N_CANDIDATES      = 25       # candidate endpoints per step

# ASF infectious period: Gamma(3, 3.5 d) → mean 10.5 d (EFSA 2018)
ASF_IP_GAMMA_SHAPE       = 3.0
ASF_IP_GAMMA_SCALE_DAYS  = 3.5
ASF_ACTIVE_STEPS_PER_DAY = 20   # active dispersal steps per day
IBMM_MAX_STEPS_CAP       = 3000 # hard per-agent cap


def simulate_ibmm(
    hab_suitability: np.ndarray,
    nodata_mask,
    source_cells: np.ndarray,
    win_tf,
    aoi_mask,
    n_agents: int = 2000,
    rng=None,
    cancel_check=None,
) -> np.ndarray:
    """Simulate n_agents infected wild-boar trajectories; return contamination-
    probability raster (visits / n_agents per cell, float32).

    Each agent starts at a random cell in source_cells, lives for a
    stochastic number of steps drawn from the ASF infectious-period
    distribution, and selects the next cell from IBMM_N_CANDIDATES
    candidates proportionally to hab_suitability at each endpoint.
    """
    if rng is None:
        rng = np.random.default_rng()

    rows, cols = hab_suitability.shape
    cellsize_m  = float(abs(win_tf.a))
    sl_scale_px = IBMM_SL_GAMMA_SCALE_M / cellsize_m

    hs = hab_suitability.astype(np.float64, copy=True)
    hs[~np.isfinite(hs)] = 0.0
    hs[hs < 0.0] = 0.0
    if nodata_mask is not None:
        hs[nodata_mask] = 0.0

    visit = np.zeros((rows, cols), dtype=np.float64)
    orr, occ = np.unravel_index(source_cells, (rows, cols))
    n_starts = len(orr)

    for agent_idx in range(n_agents):
        if cancel_check is not None and agent_idx % 100 == 0:
            if cancel_check():
                break

        life_days   = rng.gamma(ASF_IP_GAMMA_SHAPE, ASF_IP_GAMMA_SCALE_DAYS)
        agent_steps = min(
            int(round(life_days * ASF_ACTIVE_STEPS_PER_DAY)),
            IBMM_MAX_STEPS_CAP,
        )
        if agent_steps < 1:
            agent_steps = 1

        k0 = int(rng.integers(0, n_starts))
        r, c = int(orr[k0]), int(occ[k0])
        heading = float(rng.uniform(0.0, 2.0 * np.pi))

        for _ in range(agent_steps):
            if aoi_mask is not None and not aoi_mask[r, c]:
                break

            visit[r, c] += 1.0

            sls  = rng.gamma(IBMM_SL_GAMMA_SHAPE, sl_scale_px,
                             size=IBMM_N_CANDIDATES)
            tas  = rng.vonmises(0.0, IBMM_KAPPA, size=IBMM_N_CANDIDATES)
            dirs = heading + tas

            nc   = c + sls * np.cos(dirs)
            nr   = r - sls * np.sin(dirs)
            nc_i = np.round(nc).astype(np.intp)
            nr_i = np.round(nr).astype(np.intp)

            valid = ((nr_i >= 0) & (nr_i < rows) &
                     (nc_i >= 0) & (nc_i < cols))
            if not valid.any():
                break

            nr_safe = np.where(valid, nr_i, 0).astype(np.intp)
            nc_safe = np.where(valid, nc_i, 0).astype(np.intp)
            w = np.where(valid, hs[nr_safe, nc_safe], 0.0)

            wsum = float(w.sum())
            if wsum <= 0.0:
                break

            w /= wsum
            k = int(rng.choice(IBMM_N_CANDIDATES, p=w))
            if not valid[k]:
                break

            heading = float(np.arctan2(
                -(float(nr_i[k]) - r),
                float(nc_i[k]) - c,
            ))
            r, c = int(nr_i[k]), int(nc_i[k])

    prob = visit / max(n_agents, 1)
    return prob.astype(np.float32)


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

        # Outputs (each set to the GeoTIFF path when the corresponding
        # analysis step ran, None otherwise; finished() reads these to
        # decide which layers to add to the QGIS panel).
        self.risk_raster_path = None
        self.current_raster_path = None
        self.walk_raster_path = None
        self.ibmm_raster_path = None
        self.merged_lcp_polylines = []   # list[{"points", "traffic"}]
        self.forest_anchors = []         # list[((x, y), size)] - LCP targets
        self.n_lcps = 0                  # number of synthetic LCPs computed
        self.exception = None

        # Internal state
        self.nodata_mask = None
        self.aoi_mask = None
        self._cost_grid = None
        # Median resistance in the analysis window — used by the risk
        # kernel to convert ASF_DISPERSAL_KM into a matching cost scale.
        self.r_typ_for_risk = None

    # =================================================================
    # Worker thread
    # =================================================================
    def run(self):
        try:
            arr, win_tf, crs, rows, cols = self._load_window()

            # Bake fences / overpasses into the working resistance grid.
            # Plugin fences become NaN walls; plugin overpasses become
            # the minimum resistance value. Both OVERWRITE the underlying
            # raster, including pre-existing NaN walls from the notebook
            # (so an overpass on a highway makes that pixel passable).
            if self.fences or self.overpasses:
                arr = apply_modifications(
                    arr, win_tf, self.fences, self.overpasses,
                    fence_width_cells=2,
                    overpass_radius_cells=2,
                )
                # Refresh the NoData mask to reflect the post-modification
                # state. Fences turn passable cells into NaN walls; an
                # overpass turns a previously NaN cell into a finite cell.
                # Every downstream mask and percentile uses this updated
                # nodata_mask.
                self.nodata_mask = ~np.isfinite(arr)
                QgsMessageLog.logMessage(
                    f"Applied {len(self.fences)} fence(s), "
                    f"{len(self.overpasses)} overpass(es). "
                    f"Wall cells (NaN) after modifications: "
                    f"{int(self.nodata_mask.sum())}.",
                    LOG_TAG, Qgis.Info)

            # Typical resistance value for the risk-kernel scale.
            valid = ~self.nodata_mask if self.nodata_mask is not None \
                    else np.isfinite(arr)
            finite_vals = arr[valid]
            self.r_typ_for_risk = (float(np.median(finite_vals))
                                   if finite_vals.size else 1.0)
            if self.r_typ_for_risk <= 0 \
                    or not np.isfinite(self.r_typ_for_risk):
                self.r_typ_for_risk = 1.0
            QgsMessageLog.logMessage(
                f"Resistance stats: R_typical = {self.r_typ_for_risk:.2f}, "
                f"R_max = {float(np.max(finite_vals)) if finite_vals.size else 0:.2f}, "
                f"impassable cells (NaN walls): {int((~valid).sum())}.",
                LOG_TAG, Qgis.Info)

            # ---- Source disc at the outbreak pixel -------------------
            # Snaps to nearest passable cell when click lands on a wall.
            source_cells = self._point_source_cells(arr.shape, win_tf, arr=arr)
            if source_cells.size == 0:
                raise RuntimeError(
                    "No passable cells within ~2.5 km of the outbreak click. "
                    "Pick a point on land that is not entirely surrounded by "
                    "highways, large lakes, or fences.")
            QgsMessageLog.logMessage(
                f"Source disc: {source_cells.size} cell(s) at "
                f"({self.origin_xy[0]:.1f}, {self.origin_xy[1]:.1f}); "
                f"window {rows}x{cols}.",
                LOG_TAG, Qgis.Info)

            self.setProgress(10)
            if self.isCanceled():
                return False

            # ---- Single Dijkstra from outbreak point ----------------
            cost_grid, mcp = self._single_source_dijkstra(arr, source_cells)
            self._cost_grid = cost_grid
            self.setProgress(25)
            if self.isCanceled():
                return False

            # ---- Auto AOI from the cost grid -------------------------
            self.aoi_mask = self._auto_aoi_mask(cost_grid, source_cells,
                                                arr.shape)
            QgsMessageLog.logMessage(
                f"AOI: {int(self.aoi_mask.sum())} cells "
                f"({100 * self.aoi_mask.mean():.1f} % of window).",
                LOG_TAG, Qgis.Info)

            # ---- LCPs to the NEAREST low-resistance clusters ---------
            # Each LCP carries its destination cluster's size as its
            # strength; shared segments accumulate that strength so the
            # backbone corridor stands out automatically.
            if self.options.get("lcp", True):
                self._build_lcps_to_forests(
                    arr, win_tf, mcp, cost_grid,
                    n_max_targets=int(self.options.get("n_lcp_targets", 100)))
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
                        used_jl = self._try_circuitscape_jl(arr, win_tf)
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
                    self._single_source_current_scipy(arr, win_tf, crs)
            self.setProgress(80)
            if self.isCanceled():
                return False

            # ---- Biased Correlated Random Walk (optional) -----------
            if self.options.get("random_walk", False):
                self._random_walks(
                    arr, win_tf, crs, source_cells,
                    n_walks=int(self.options.get("n_walks", 500)),
                    kappa=float(self.options.get("walk_kappa", 2.0)),
                )
            self.setProgress(90)
            if self.isCanceled():
                return False

            # ---- iSSF-IBMM agent simulation (optional) --------------
            if self.options.get("ibmm", False):
                self._run_ibmm(arr, win_tf, source_cells)
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

    def _snap_to_passable(self, arr, r0, c0, max_radius_cells=100):
        """Find the (r, c) of the nearest passable cell to (r0, c0).

        Searches a square window of side 2*max_radius_cells (default
        ~2.5 km at 25 m/px) and returns the Euclidean-nearest cell with
        a finite resistance. Returns None if no such cell exists in the
        search window.
        """
        rows, cols = arr.shape
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
        dr = rr - (rc - r1)
        dc = ccs - (cc - c1)
        d2 = dr * dr + dc * dc
        k = int(np.argmin(d2))
        return (int(rr[k] + r1), int(ccs[k] + c1))

    # -----------------------------------------------------------------
    def _single_source_dijkstra(self, arr, source_cells):
        """One Dijkstra sweep from the source disc through the whole grid.

        skimage.graph.MCP_Geometric crashes on NaN; substitute NaN → inf
        for this call only. The original arr (with NaN walls) is kept for
        everything else.
        """
        from skimage.graph import MCP_Geometric
        rows, cols = arr.shape
        rr, cc = np.unravel_index(source_cells, (rows, cols))
        starts = list(zip(rr.tolist(), cc.tolist()))

        arr_for_mcp = np.where(np.isfinite(arr), arr, np.inf)
        mcp = MCP_Geometric(arr_for_mcp, fully_connected=True)
        cost_grid, _ = mcp.find_costs(starts)
        return cost_grid, mcp

    # -----------------------------------------------------------------
    def _auto_aoi_mask(self, cost_grid, source_cells, shape):
        """Every cell reachable from the source disc.

        AOI = cells with finite Dijkstra cost. Cells behind NaN walls
        get cost == inf and are excluded automatically.
        """
        rows, cols = shape
        aoi = np.isfinite(cost_grid)
        rr, cc = np.unravel_index(source_cells, (rows, cols))
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
    # Rasterise the user-supplied habitat polygons onto the analysis
    # grid. Each polygon becomes one labelled cluster (label 1, 2, ...),
    # so the rest of the LCP step can treat user habitats and
    # auto-detected clusters identically.
    # -----------------------------------------------------------------
    def _rasterize_user_habitats(self, habitats_geojson, win_tf, shape):
        if not habitats_geojson:
            return None
        rows, cols = shape
        labeled = np.zeros((rows, cols), dtype=np.int32)
        # Rasterise each polygon individually so we get one label
        # per habitat (rasterio.features.rasterize on a list with
        # different values does this in one call).
        shapes = [(g, i + 1) for i, g in enumerate(habitats_geojson)]
        labeled = rasterize(
            shapes,
            out_shape=shape,
            transform=win_tf,
            fill=0,
            dtype="int32",
            all_touched=True,
        )
        return labeled

    # -----------------------------------------------------------------
    # LCPs to the NEAREST low-resistance clusters, weighted by cluster
    # size.
    #
    # Method:
    #   1. Threshold the in-AOI valid cells at the 20th percentile to
    #      isolate "low-resistance" pixels (good boar habitat).
    #   2. Find connected components (scipy.ndimage.label, 8-connected).
    #      Cluster SIZE (cell count) is its strength / importance:
    #      more habitat there means more boars potentially dispersing
    #      to and from it.
    #   3. Skip the cluster(s) overlapping the origin disc.
    #   4. For each remaining cluster, pick the cell with the LOWEST
    #      Dijkstra cost from the origin as the anchor (the natural
    #      entry point), and record (entry_cost, anchor_r, anchor_c,
    #      cluster_size).
    #   5. Sort by cost ASCENDING (nearest in boar-cost terms first),
    #      cap at n_max_targets (default 100).
    #   6. Traceback each anchor through the MCP.
    #   7. Edge strength: each LCP carries its destination cluster's
    #      SIZE. For every edge along the path, ADD that size to the
    #      edge's cumulative strength. Where multiple LCPs share a
    #      segment, the segment's strength is the SUM of every
    #      destination it serves - the shared backbone matters more
    #      because more habitat depends on it.
    #   8. Merge consecutive edges with the same cumulative strength
    #      into single line features.
    # -----------------------------------------------------------------
    def _build_lcps_to_forests(self, arr, win_tf, mcp, cost_grid,
                               n_max_targets=100,
                               forest_quantile=0.20,
                               min_patch_cells=5):
        rows, cols = arr.shape

        # --- Choose the destination source ---------------------------
        # Path A: user-supplied core-habitat shapefile.
        # Path B: auto-detected low-resistance clusters (fallback when
        # no habitat layer is selected in the dialog).
        user_habitats = self.options.get("habitats_geojson")
        if user_habitats:
            labeled = self._rasterize_user_habitats(
                user_habitats, win_tf, arr.shape)
            if labeled is None or int(labeled.max()) == 0:
                QgsMessageLog.logMessage(
                    "Habitat polygons rasterise to zero cells; "
                    "LCPs skipped.",
                    LOG_TAG, Qgis.Warning)
                return
            n_features = int(labeled.max())
            QgsMessageLog.logMessage(
                f"LCP target source: user-supplied habitat polygons "
                f"({n_features} habitats rasterised on the grid).",
                LOG_TAG, Qgis.Info)
        else:
            try:
                from scipy.ndimage import label
            except ImportError:
                QgsMessageLog.logMessage(
                    "scipy.ndimage missing; LCPs skipped.",
                    LOG_TAG, Qgis.Warning)
                return
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
            labeled, n_features = label(
                forest_mask, structure=np.ones((3, 3), dtype=int))
            if n_features == 0:
                return
            QgsMessageLog.logMessage(
                f"LCP target source: auto-detected low-resistance "
                f"clusters ({n_features} clusters).",
                LOG_TAG, Qgis.Info)

        # --- Identify the cluster(s) overlapping the source disc -----
        origin_ids = self._point_source_cells((rows, cols), win_tf)
        origin_patch_ids = set()
        if origin_ids.size:
            orr, occ = np.unravel_index(origin_ids, (rows, cols))
            lab_at_origin = labeled[orr, occ]
            origin_patch_ids = {int(v) for v in lab_at_origin if v > 0}

        # --- Anchor cell per (non-origin) cluster --------------------
        # Anchor = the cell in the cluster with the LOWEST Dijkstra
        # cost from the origin. Conceptually this is the cluster's
        # "natural entry point" for a dispersing boar, not its
        # geometric centroid.
        clusters = []   # list of (entry_cost, anchor_r, anchor_c, size)
        for cid in range(1, n_features + 1):
            if cid in origin_patch_ids:
                continue
            rr_p, cc_p = np.where(labeled == cid)
            size = rr_p.size
            if size < min_patch_cells:
                continue
            costs = cost_grid[rr_p, cc_p]
            finite = np.isfinite(costs)
            if not finite.any():
                continue
            k = int(np.argmin(np.where(finite, costs, np.inf)))
            entry_cost = float(costs[k])
            if not np.isfinite(entry_cost):
                continue
            clusters.append((entry_cost,
                             int(rr_p[k]), int(cc_p[k]), int(size)))

        if not clusters:
            QgsMessageLog.logMessage(
                "Only the origin's own cluster is in the AOI; no LCPs.",
                LOG_TAG, Qgis.Info)
            return

        # Nearest clusters first (by boar cost-distance), cap at n_max_targets.
        clusters.sort(key=lambda x: x[0])
        clusters = clusters[:int(n_max_targets)]
        QgsMessageLog.logMessage(
            f"LCP destinations: {len(clusters)} nearest low-resistance "
            f"clusters. First 5 sizes: {[c[3] for c in clusters[:5]]}. "
            f"First 5 entry costs: {[round(c[0], 1) for c in clusters[:5]]}.",
            LOG_TAG, Qgis.Info)

        # Remember anchors for the visualisation layer (in world coords).
        self.forest_anchors = [
            (rio_xy(win_tf, r, c), size)
            for (_cost, r, c, size) in clusters
        ]

        # --- Traceback each anchor and carry its cluster size --------
        paths_with_strength = []   # list of (cells, cluster_size)
        for _cost, r, c, size in clusters:
            try:
                indices = mcp.traceback((r, c))
            except Exception:
                continue
            cells = [(int(ri), int(ci)) for ri, ci in indices]
            if len(cells) >= 2:
                paths_with_strength.append((cells, size))
        self.n_lcps = len(paths_with_strength)
        if not paths_with_strength:
            return

        # --- Cumulative edge strength --------------------------------
        # Each edge's strength = SUM of destination cluster sizes of
        # every LCP that traverses it. Shared backbone segments stack
        # up automatically: a corridor leading to several large
        # clusters ends up far heavier than a one-cluster detour.
        from collections import defaultdict
        edge_strength = defaultdict(int)
        for cells, strength in paths_with_strength:
            for a, b in zip(cells, cells[1:]):
                edge = (a, b) if a < b else (b, a)
                edge_strength[edge] += int(strength)

        # --- Merge into uniform-strength polylines -------------------
        visited = set()
        runs = []
        for cells, _strength in paths_with_strength:
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
                s = edge_strength[edge]
                if seg_count is None:
                    seg_cells = [a, b]
                    seg_count = s
                elif s == seg_count:
                    seg_cells.append(b)
                else:
                    runs.append((seg_cells, seg_count))
                    seg_cells = [a, b]
                    seg_count = s
            if seg_count is not None and len(seg_cells) >= 2:
                runs.append((seg_cells, seg_count))

        for cells, strength in runs:
            pts = [QgsPointXY(*rio_xy(win_tf, r, c)) for r, c in cells]
            self.merged_lcp_polylines.append({
                "points":  pts,
                "traffic": int(strength),   # field name kept; value =
                                            # cumulative cluster strength
            })

    # =================================================================
    # Output rasters
    # =================================================================
    def _build_risk_raster(self, cost_grid, win_tf):
        """Continuous infection-risk raster: exp(-cost / D_cost).

        The decay scale D_cost is the cost of traversing ASF_DISPERSAL_KM
        through typical-resistance terrain. This pins the kernel to the
        published mean wild boar dispersal distance (4 km, see citations
        at the top of this module) instead of letting the AOI size
        stretch the bands outward.
        """
        cellsize_m = float(abs(win_tf.a))
        r_typ = float(self.r_typ_for_risk or 1.0)
        # Number of pixels to walk to cover ASF_DISPERSAL_KM, and the
        # approximate cumulative cost along that walk through typical
        # terrain.
        n_pixels_D = ASF_DISPERSAL_KM * 1000.0 / cellsize_m
        scale = max(n_pixels_D * r_typ, 1e-9)

        QgsMessageLog.logMessage(
            f"Risk kernel: D = {ASF_DISPERSAL_KM:.1f} km, cellsize = "
            f"{cellsize_m:.0f} m, R_typ = {r_typ:.2f}, cost-scale = "
            f"{scale:.1f}. risk(d) = exp(-cost / cost-scale).",
            LOG_TAG, Qgis.Info)

        risk = np.exp(-cost_grid / scale)
        risk[~np.isfinite(cost_grid)] = np.nan
        return self._mask_crop_write_raster(
            risk.astype(np.float32), win_tf, prefix="risk")

    # -----------------------------------------------------------------
    def _point_source_cells(self, shape, win_tf, arr=None):
        """Flat-index array of a small disc (radius = origin_radius_cells)
        around self.origin_xy.

        When arr is supplied, cells that are impassable (NaN) are stripped
        from the disc. If the entire disc is impassable, the centre snaps
        to the nearest passable cell via _snap_to_passable and self.origin_xy
        is updated accordingly.
        """
        rows, cols = shape
        try:
            inv = ~win_tf
        except Exception:
            return np.array([], dtype=np.intp)
        c0, r0 = inv * self.origin_xy
        r0_i = int(round(r0))
        c0_i = int(round(c0))
        radius = max(1, self.origin_radius_cells)

        rr, cc = np.meshgrid(
            np.arange(max(0, r0_i - radius), min(rows, r0_i + radius + 1)),
            np.arange(max(0, c0_i - radius), min(cols, c0_i + radius + 1)),
            indexing="ij",
        )
        in_disc = (rr - r0_i) ** 2 + (cc - c0_i) ** 2 <= radius * radius
        rr_d = rr[in_disc].astype(np.intp)
        cc_d = cc[in_disc].astype(np.intp)

        if arr is not None and rr_d.size > 0:
            passable = np.isfinite(arr[rr_d, cc_d])
            if passable.any():
                rr_d = rr_d[passable]
                cc_d = cc_d[passable]
            else:
                snap = self._snap_to_passable(arr, r0_i, c0_i)
                if snap is None:
                    return np.array([], dtype=np.intp)
                r0_i, c0_i = snap
                x_snap, y_snap = rio_xy(win_tf, r0_i, c0_i)
                QgsMessageLog.logMessage(
                    "Origin click is on an impassable cell. Snapped to "
                    f"nearest passable cell: "
                    f"({self.origin_xy[0]:.1f}, {self.origin_xy[1]:.1f})"
                    f" → ({x_snap:.1f}, {y_snap:.1f}).",
                    LOG_TAG, Qgis.Warning)
                self.origin_xy = (float(x_snap), float(y_snap))
                return self._point_source_cells(shape, win_tf, arr=arr)

        if rr_d.size == 0:
            if 0 <= r0_i < rows and 0 <= c0_i < cols:
                return np.array([r0_i * cols + c0_i], dtype=np.intp)
            return np.array([], dtype=np.intp)
        return (rr_d * cols + cc_d).astype(np.intp)

    # Pinchpoint raster: Circuit theory with AOI boundary as sink.
    # -----------------------------------------------------------------
    def _try_circuitscape_jl(self, arr, win_tf):
        rows, cols = arr.shape
        cs_source_ids = self._point_source_cells(arr.shape, win_tf)
        if cs_source_ids.size == 0:
            return False
        source_mask = np.zeros((rows, cols), dtype=bool)
        srr, scc = np.unravel_index(cs_source_ids, (rows, cols))
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

        # Suppress the source zone: all source cells are held at the
        # injection voltage, so the voltage gradient between them is
        # nearly zero → near-zero current density → misleading dark
        # "dead zone" at the outbreak origin on the Magma ramp.
        # The corridor map should show where disease spreads TO, not
        # that the origin exists (that's already shown by the origin
        # marker). Setting source cells to NaN removes the artifact
        # without affecting any corridor information.
        current = current.astype(np.float32)
        current[source_mask] = np.nan

        self.current_raster_path = self._mask_crop_write_raster(
            current, win_tf, prefix="pinchpoints_cs")
        return True

    def _single_source_current_scipy(self, arr, win_tf, crs):
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

        cs_source_ids = self._point_source_cells(arr.shape, win_tf)
        if cs_source_ids.size == 0:
            return
        bvec = np.zeros(n)
        bvec[cs_source_ids] += 1.0 / cs_source_ids.size

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
        node_current = node_current.reshape(rows, cols).astype(np.float32)

        # Suppress the point-source disc (same NaN masking as the
        # Circuitscape.jl path).
        srr, scc = np.unravel_index(cs_source_ids, (rows, cols))
        node_current[srr, scc] = np.nan

        self.current_raster_path = self._mask_crop_write_raster(
            node_current, win_tf, prefix="pinchpoints")

    # -----------------------------------------------------------------
    # Biased Correlated Random Walk (BCRW) from the origin.
    #
    # Each step samples one of the 8 neighbours with probability
    #
    #     P(neighbour) proportional to   w_cost * w_dir
    #
    # where the two factors capture the two well-established
    # behavioural components of large-mammal movement:
    #
    #     w_cost = 1 / (R_neighbour * step_distance)
    #
    #         Boars prefer easier terrain and avoid wasting energy
    #         on longer (diagonal) steps. R_neighbour is the
    #         post-penalty resistance; step_distance is 1 for
    #         orthogonal moves and sqrt(2) for diagonals.
    #
    #     w_dir  = exp(kappa * cos(theta_step - theta_heading))
    #
    #         A von Mises directional-persistence term. After a
    #         step in direction theta_heading the walker is more
    #         likely to keep going in roughly the same direction.
    #         kappa = 0 collapses to an unbiased random walk;
    #         kappa -> infinity collapses to deterministic
    #         straight-line motion. We use kappa = 2 (moderate
    #         persistence; turning-angle SD ~ 50 deg), consistent
    #         with the cos_ta coefficient from the in-matrix iSSF
    #         in Step 5 of the resistance-surface notebook.
    #
    # References:
    #     Codling, E. A., Plank, M. J., Benhamou, S. (2008).
    #         Random walk models in biology. Journal of the Royal
    #         Society Interface 5(25): 813-834.
    #     Turchin, P. (1998). Quantitative Analysis of Movement.
    #         Sinauer.
    #     Avgar, T., Potts, J. R., Lewis, M. A., Boyce, M. S.
    #         (2016). Integrated step selection analysis. Methods
    #         in Ecology and Evolution 7(5): 619-630.
    #
    # The output raster is the cell-visit count over n_walks walks.
    # In the limit n_walks -> infinity it converges to the same
    # surface as the Circuit-theory current map (Saerens et al.
    # 2009, RSP framework), which is a useful sanity check.
    # -----------------------------------------------------------------
    def _random_walks(self, arr, win_tf, crs, source_cells,
                      n_walks=500, max_steps=2000, kappa=2.0):
        rows, cols = arr.shape
        rr, cc = np.unravel_index(source_cells, (rows, cols))
        starts = np.column_stack([rr, cc])

        visit = np.zeros((rows, cols), dtype=np.int64)
        rng = np.random.default_rng()

        # 8-neighbour table: (dr, dc), step distance, step direction
        # as a unit vector in (dx, dy) raster coords.
        offsets = np.array([
            (-1, -1), (-1, 0), (-1, 1),
            ( 0, -1),          ( 0, 1),
            ( 1, -1), ( 1, 0), ( 1, 1),
        ], dtype=np.intp)
        SQ2 = float(np.sqrt(2.0))
        step_dist = np.array([
            SQ2, 1.0, SQ2,
            1.0,      1.0,
            SQ2, 1.0, SQ2,
        ], dtype=np.float64)
        # Unit vector pointing from current cell to neighbour cell.
        # dx is the column delta, dy the row delta (raster row indexes
        # increase downward, which does not matter for cosine of angle
        # differences).
        ux = np.array([-1, 0, 1, -1, 1, -1, 0, 1], dtype=np.float64)
        uy = np.array([-1,-1,-1,  0, 0,  1, 1, 1], dtype=np.float64)
        norm = np.sqrt(ux * ux + uy * uy)
        ux /= norm; uy /= norm

        aoi    = self.aoi_mask
        nodata = self.nodata_mask

        for _ in range(int(n_walks)):
            if self.isCanceled():
                return

            # Random start inside the origin disc, random initial heading
            # (uniform over [0, 2 pi]).
            r, c = starts[rng.integers(0, len(starts))]
            r, c = int(r), int(c)
            theta = float(rng.uniform(0.0, 2.0 * np.pi))
            hx, hy = float(np.cos(theta)), float(np.sin(theta))

            for _ in range(int(max_steps)):
                visit[r, c] += 1

                # Stop when the walker leaves the AOI.
                if aoi is not None and not aoi[r, c]:
                    break

                nrs = r + offsets[:, 0]
                ncs = c + offsets[:, 1]
                ok = (nrs >= 0) & (nrs < rows) \
                   & (ncs >= 0) & (ncs < cols)

                # NoData / walls are forbidden moves.
                if nodata is not None and ok.any():
                    safe_r = np.clip(nrs, 0, rows - 1)
                    safe_c = np.clip(ncs, 0, cols - 1)
                    ok = ok & ~nodata[safe_r, safe_c]

                if not ok.any():
                    break

                nrs_ok = nrs[ok]
                ncs_ok = ncs[ok]
                d_ok   = step_dist[ok]
                ux_ok  = ux[ok]
                uy_ok  = uy[ok]

                # Cost bias.
                w_cost = 1.0 / (arr[nrs_ok, ncs_ok] * d_ok)

                # Directional bias: cosine of the angle between the
                # current heading (hx, hy) and the step direction
                # (ux, uy), pushed through the von Mises exponential.
                cos_angle = hx * ux_ok + hy * uy_ok
                w_dir = np.exp(float(kappa) * cos_angle)

                w = w_cost * w_dir
                w_sum = w.sum()
                if not np.isfinite(w_sum) or w_sum <= 0:
                    break
                w /= w_sum

                k = int(rng.choice(len(nrs_ok), p=w))
                r, c = int(nrs_ok[k]), int(ncs_ok[k])
                # Update the heading to the actual direction taken.
                hx, hy = float(ux_ok[k]), float(uy_ok[k])

        if visit.any():
            grid = visit.astype(np.float32)
            grid[grid == 0] = np.nan
            self.walk_raster_path = self._mask_crop_write_raster(
                grid, win_tf, prefix="random_walk")
            QgsMessageLog.logMessage(
                f"BCRW: {n_walks} walks, max {max_steps} steps each, "
                f"kappa = {kappa:.1f}. Total cell visits: "
                f"{int(visit.sum())}.",
                LOG_TAG, Qgis.Info)

    # =================================================================
    # iSSF-IBMM: agent simulation using habitat suitability as the
    # step-selection surface.
    # =================================================================
    def _run_ibmm(self, arr, win_tf, source_cells):
        """Simulate n_agents infected wild-boar trajectories (iSSF-IBMM).

        Each agent:
          - Starts at a random cell in the outbreak zone.
          - Lives for a stochastic number of active steps sampled from
            Gamma(ASF_IP_GAMMA_SHAPE, ASF_IP_GAMMA_SCALE_DAYS) × steps/day,
            representing the time from infection to death (EFSA 2018).
          - At each step selects the next cell from IBMM_N_CANDIDATES
            candidates proportional to 1/resistance (the iSSF selection rule).
        Output: contamination-probability raster = fraction of agents that
        visit each cell.
        """
        n_agents = int(self.options.get("n_agents", 2000))

        # Recover habitat suitability S from resistance R using the exact
        # inverse of the Keeley exponential transformation applied in
        # scripts/01_create_resistance_surface.ipynb:
        #
        #   Forward:  R = exp(C * (1 - S))    C = 4
        #   Inverse:  S = 1 - log(R) / C
        #
        # Range: R ∈ [1, exp(4)≈54.6]  →  S ∈ [1, 0]
        # Impassable cells (NaN / negative) are set to S = 0.
        KEELEY_C = 4.0
        finite_mask = np.isfinite(arr) & (arr > 0)
        hs = np.where(finite_mask, 1.0 - np.log(arr) / KEELEY_C, 0.0)
        hs = np.clip(hs, 0.0, 1.0).astype(np.float64)

        mean_ip_days = ASF_IP_GAMMA_SHAPE * ASF_IP_GAMMA_SCALE_DAYS
        mean_steps = int(round(mean_ip_days * ASF_ACTIVE_STEPS_PER_DAY))
        QgsMessageLog.logMessage(
            f"iSSF-IBMM: {n_agents} agents; "
            f"infectious period Gamma(shape={ASF_IP_GAMMA_SHAPE}, "
            f"scale={ASF_IP_GAMMA_SCALE_DAYS} d) "
            f"→ mean {mean_ip_days:.1f} d / {mean_steps} steps; "
            f"step Gamma(shape={IBMM_SL_GAMMA_SHAPE}, "
            f"scale={IBMM_SL_GAMMA_SCALE_M} m); "
            f"kappa={IBMM_KAPPA}; "
            f"{IBMM_N_CANDIDATES} candidates/step.",
            LOG_TAG, Qgis.Info)

        prob = simulate_ibmm(
            hab_suitability=hs,
            nodata_mask=self.nodata_mask,
            source_cells=source_cells,
            win_tf=win_tf,
            aoi_mask=self.aoi_mask,
            n_agents=n_agents,
            cancel_check=self.isCanceled,
        )

        if prob is not None and float(prob.sum()) > 0:
            grid = prob.copy()

            # Threshold: drop cells visited by fewer than 1 agent.
            # This scales with n_agents (principled minimum detectable signal)
            # and avoids showing numerical noise from near-zero visitation.
            min_prob = 1.0 / max(n_agents, 1)
            grid[grid < min_prob] = np.nan

            # Log₁₀ transform: contamination probability spans ~3 orders of
            # magnitude (origin ≈ 1, far cells ≈ 1/n_agents ≈ 5e-4). Linear
            # scale collapses all distant structure into near-black; log scale
            # spreads the gradient across the full color ramp so corridors and
            # barriers are readable everywhere, not just at the outbreak zone.
            valid = np.isfinite(grid)
            if valid.any():
                grid[valid] = np.log10(grid[valid])
            # Values now range from log10(1/n_agents) ≈ -3.3 to 0.0

            self.ibmm_raster_path = self._mask_crop_write_raster(
                grid, win_tf, prefix="ibmm")
            peak_raw = float(np.nanmax(prob)) if prob is not None else 0.0
            QgsMessageLog.logMessage(
                f"iSSF-IBMM: completed. Peak contamination probability: "
                f"{peak_raw:.4f} ({peak_raw * 100:.1f}% of agents). "
                f"Output: log10(p), range "
                f"[{float(np.nanmin(grid)):.2f}, {float(np.nanmax(grid)):.2f}] "
                f"(threshold: p ≥ {min_prob:.4f} = 1 agent).",
                LOG_TAG, Qgis.Info)

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
                               "Dispersal corridor probability (Circuit theory, raw current density)")
            if r.isValid():
                self._apply_singleband_ramp(r, "Magma")
                add_wildboar_layer(r, ZOrder.PINCHPOINTS)

        # Continuous infection-risk raster
        if self.risk_raster_path:
            r = QgsRasterLayer(self.risk_raster_path,
                               "Infection risk (exponential decay function)")
            if r.isValid():
                self._apply_risk_ramp(r)
                add_wildboar_layer(r, ZOrder.RISK)

        # Biased Correlated Random Walk density
        if self.walk_raster_path:
            n = int(self.options.get("n_walks", 500))
            r = QgsRasterLayer(
                self.walk_raster_path,
                f"Random walk density (BCRW, n = {n})")
            if r.isValid():
                self._apply_singleband_ramp(r, "Viridis")
                add_wildboar_layer(r, ZOrder.WALK_DENSITY)

        # iSSF-IBMM contamination probability (log10-scaled)
        if self.ibmm_raster_path:
            n = int(self.options.get("n_agents", 2000))
            r = QgsRasterLayer(
                self.ibmm_raster_path,
                f"ASF contamination probability log₁₀(p) "
                f"— iSSF-IBMM, n = {n} agents")
            if r.isValid():
                self._apply_singleband_ramp(r, "Magma")
                add_wildboar_layer(r, ZOrder.IBMM_DENSITY)

        QgsMessageLog.logMessage(
            "Done. lcps={l}({s} segments) pinch={p} risk={r} walks={w}".format(
                l=self.n_lcps,
                s=len(self.merged_lcp_polylines),
                p="yes" if self.current_raster_path else "no",
                r="yes" if self.risk_raster_path else "no",
                w="yes" if self.walk_raster_path else "no"),
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
            f"Dispersal corridors (LCPs to {self.n_lcps} nearest "
            f"clusters, weighted by cluster size; {len(polys)} segments)",
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
                             f"Weak (strength < {int(low_hi)})"),
            QgsRendererRange(low_hi, med_hi,
                             line_sym("#2874A6", 1.6),
                             f"Medium ({int(low_hi)}-{int(med_hi)})"),
            QgsRendererRange(med_hi, upper,
                             line_sym("#154360", 3.2),
                             f"Strong backbone (>= {int(med_hi)})"),
        ]
        renderer = QgsGraduatedSymbolRenderer("traffic", ranges)
        layer.setRenderer(renderer)
        add_wildboar_layer(layer, ZOrder.LCP_TRAFFIC)

    # -----------------------------------------------------------------
    def _add_forest_anchors_layer(self):
        """Show the auto-detected forest destinations as labelled dots."""
        layer = QgsVectorLayer(
            f"Point?crs={self.crs_wkt}",
            f"Cluster destinations ({len(self.forest_anchors)} "
            f"nearest low-resistance clusters)",
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
        add_wildboar_layer(layer, ZOrder.ANCHORS)

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
