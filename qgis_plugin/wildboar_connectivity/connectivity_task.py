# -*- coding: utf-8 -*-
"""
Background QgsTask that computes:

    * Pairwise least-cost paths (Dijkstra via scikit-image) between every
      pair of selected core habitats.

    * A cumulative Circuitscape-style "current flow" raster summing the
      absolute edge currents from every pairwise solve. Bright pixels in
      this map are pinchpoints: cells where many corridors converge and
      whose loss would disconnect the network.

Math summary (full derivation in main_plugin.py docstring):
    g_ij = 1 / (0.5 * (R_i + R_j))            # conductance
    L    = D - A    (weighted graph Laplacian, D=sum of incident g)
    L v = b   with  b[source_cells] = +1/n_s, b[sink_cells] = -1/n_s.
    I_e  = g_ij * (v_i - v_j)                 # edge current
    Cum_e = sum over all habitat pairs of |I_e|.
    Node value = 0.5 * sum_{j ~ i} Cum(i,j).

Performance: the Laplacian is factored ONCE (scipy splu) with a single
corner cell grounded; every pair is then a fast back-solve.
"""

import os
import tempfile
import itertools
import json

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.windows import from_bounds, Window
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


LOG_TAG = "wildboar-v2"


class HabitatConnectivityTask(QgsTask):
    """Compute pairwise LCPs and a cumulative current-flow map across all
    habitats falling inside a user-defined wild boar range buffer.

    Constructor arguments are plain Python data (lists, dicts, strings),
    deliberately not QGIS objects, because QGIS layers are NOT thread-safe.
    The UI thread extracts everything before this task is launched.
    """

    def __init__(self,
                 raster_path: str,
                 habitats: list,         # [{id, label, geojson, area}]
                 target_id: int,
                 window_bounds: tuple,   # (min_x, min_y, max_x, max_y) raster CRS
                 crs_wkt: str,           # raster CRS for output layers
                 run_lcp: bool,
                 run_circuit: bool,
                 ramp_name: str = "Magma"):
        super().__init__("Wildboar habitat connectivity", QgsTask.CanCancel)
        self.raster_path = raster_path
        self.habitats = habitats
        self.target_id = target_id
        self.window_bounds = window_bounds
        self.crs_wkt = crs_wkt
        self.run_lcp = run_lcp
        self.run_circuit = run_circuit
        self.ramp_name = ramp_name

        # Outputs (filled by run(), read by finished())
        self.lcp_records = []          # list of {from_id, to_id, cost, points}
        self.current_raster_path = None
        self.n_pairs = 0
        self.exception = None

    # =================================================================
    # Worker thread
    # =================================================================
    def run(self):
        try:
            arr, win_tf, crs, rows, cols, nodata_mask = self._load_window()

            # Rasterize every habitat polygon onto the analysis grid.
            # `flat_ids` is the 1-D index list for each habitat.
            habitat_masks = self._rasterize_habitats(arr.shape, win_tf)
            if len(habitat_masks) < 2:
                raise RuntimeError(
                    f"Only {len(habitat_masks)} habitats rasterise onto the "
                    f"raster window. Need at least 2.")

            pairs = list(itertools.combinations(habitat_masks.keys(), 2))
            self.n_pairs = len(pairs)
            QgsMessageLog.logMessage(
                f"Analysing {len(habitat_masks)} habitats, {self.n_pairs} pairs "
                f"on a {rows}x{cols} window.", LOG_TAG, Qgis.Info)

            if self.isCanceled():
                return False

            # ---- LCPs (always cheap; skip only if disabled) -----------
            if self.run_lcp:
                self._compute_lcps(arr, win_tf, habitat_masks, pairs)
            self.setProgress(35 if self.run_circuit else 95)

            if self.isCanceled():
                return False

            # ---- Cumulative current flow ------------------------------
            if self.run_circuit:
                self._compute_cumulative_currents(
                    arr, win_tf, crs, habitat_masks, pairs)
            self.setProgress(100)
            return True

        except Exception as exc:
            self.exception = exc
            QgsMessageLog.logMessage(f"Task failed: {exc}",
                                     LOG_TAG, Qgis.Critical)
            return False

    # -----------------------------------------------------------------
    def _load_window(self):
        """Read the resistance raster within the analysis window.

        NoData is replaced with a very large (but finite) resistance so the
        Laplacian stays connected; the NoData mask is returned so we can
        respect it later if needed.
        """
        min_x, min_y, max_x, max_y = self.window_bounds
        with rasterio.open(self.raster_path) as src:
            crs = src.crs
            win = from_bounds(min_x, min_y, max_x, max_y, src.transform)
            win = win.intersection(Window(0, 0, src.width, src.height))
            if win.width <= 0 or win.height <= 0:
                raise RuntimeError("Analysis window does not overlap raster.")
            arr = src.read(1, window=win).astype(np.float64)
            win_tf = src.window_transform(win)
            nodata = src.nodata

        nodata_mask = np.zeros(arr.shape, dtype=bool)
        if nodata is not None:
            nodata_mask = (arr == nodata)
            arr[nodata_mask] = 1e6
        # Keep resistance in a numerically sane range.
        arr = np.clip(arr, 1e-3, 1e6)
        rows, cols = arr.shape
        return arr, win_tf, crs, rows, cols, nodata_mask

    # -----------------------------------------------------------------
    def _rasterize_habitats(self, shape, win_tf):
        """Return {habitat_id: np.ndarray of flat cell indices in the window}."""
        out = {}
        for h in self.habitats:
            geom = h["geojson"]
            mask = rasterize(
                [(geom, 1)],
                out_shape=shape,
                transform=win_tf,
                fill=0,
                dtype="uint8",
                all_touched=True,
            )
            flat = np.flatnonzero(mask.flatten())
            if flat.size:
                out[h["id"]] = flat
        return out

    # -----------------------------------------------------------------
    def _habitat_anchor_rc(self, flat_ids, rows, cols):
        """Pick a representative (row, col) cell inside a habitat.

        Uses the centroid of the rasterised cells if it's actually inside
        the habitat, otherwise falls back to the nearest inside cell.
        """
        rr, cc = np.unravel_index(flat_ids, (rows, cols))
        cr = int(round(rr.mean()))
        cc_ = int(round(cc.mean()))
        if cr * cols + cc_ in set(flat_ids.tolist()):
            return cr, cc_
        d2 = (rr - cr) ** 2 + (cc - cc_) ** 2
        k = int(np.argmin(d2))
        return int(rr[k]), int(cc[k])

    # -----------------------------------------------------------------
    def _compute_lcps(self, arr, win_tf, habitat_masks, pairs):
        from skimage.graph import route_through_array
        rows, cols = arr.shape

        anchors = {hid: self._habitat_anchor_rc(idx, rows, cols)
                   for hid, idx in habitat_masks.items()}

        for k, (a, b) in enumerate(pairs):
            if self.isCanceled():
                return
            try:
                indices, cost = route_through_array(
                    arr, anchors[a], anchors[b],
                    fully_connected=True, geometric=True)
                pts = [QgsPointXY(*rio_xy(win_tf, r, c)) for r, c in indices]
                self.lcp_records.append({
                    "from_id": int(a), "to_id": int(b),
                    "cost": float(cost), "points": pts,
                })
            except Exception as exc:
                QgsMessageLog.logMessage(
                    f"LCP {a}<->{b} failed: {exc}", LOG_TAG, Qgis.Warning)
            self.setProgress(5 + 25 * (k + 1) / max(1, self.n_pairs))

    # -----------------------------------------------------------------
    def _compute_cumulative_currents(self, arr, win_tf, crs,
                                     habitat_masks, pairs):
        """Build the Laplacian once, factor once, solve once per pair.

        The trick: ground a single (constant) corner cell so the reduced
        Laplacian is non-singular and can be LU-factored. The factorisation
        depends only on the resistance surface, not the source/sink
        choice, so every pairwise solve is a cheap back-substitution.
        """
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import splu

        rows, cols = arr.shape
        n = rows * cols
        R_flat = arr.flatten()

        # Edge endpoint arrays (4-connected grid).
        idx = np.arange(n).reshape(rows, cols)
        h_i = idx[:, :-1].flatten();  h_j = idx[:, 1:].flatten()
        v_i = idx[:-1, :].flatten();  v_j = idx[1:, :].flatten()
        e_i = np.concatenate([h_i, v_i])
        e_j = np.concatenate([h_j, v_j])

        # Conductance from resistance (cells in series, half each).
        R_edge = 0.5 * (R_flat[e_i] + R_flat[e_j])
        g_edge = 1.0 / np.maximum(R_edge, 1e-9)

        # Weighted Laplacian L = D - A in CSC.
        data = np.concatenate([-g_edge, -g_edge,  g_edge,  g_edge])
        row  = np.concatenate([   e_i,    e_j,    e_i,    e_j])
        col  = np.concatenate([   e_j,    e_i,    e_i,    e_j])
        L = coo_matrix((data, (row, col)), shape=(n, n)).tocsc()

        # Ground one cell (corner) and factor ONCE.
        ground = 0
        keep = np.ones(n, dtype=bool)
        keep[ground] = False
        L_red = L[keep][:, keep]
        self.setProgress(45)
        lu = splu(L_red)
        self.setProgress(55)

        cum_edge = np.zeros(e_i.size, dtype=np.float64)

        for k, (a, b) in enumerate(pairs):
            if self.isCanceled():
                return
            a_ids = habitat_masks[a]
            b_ids = habitat_masks[b]
            bvec = np.zeros(n)
            bvec[a_ids] += 1.0 / a_ids.size       # +1 A spread over A
            bvec[b_ids] -= 1.0 / b_ids.size       # -1 A spread over B
            b_red = bvec[keep]

            v_red = lu.solve(b_red)
            v = np.zeros(n)
            v[keep] = v_red

            edge_I = g_edge * (v[e_i] - v[e_j])
            cum_edge += np.abs(edge_I)

            self.setProgress(55 + 40 * (k + 1) / max(1, self.n_pairs))

        # Node current = 0.5 * sum of |edge currents| incident on node.
        node_current = np.zeros(n)
        np.add.at(node_current, e_i, cum_edge)
        np.add.at(node_current, e_j, cum_edge)
        node_current *= 0.5
        node_current = node_current.reshape(rows, cols)

        # Mask habitat cells: they are core habitat, not corridor, and the
        # injection geometry biases their values upward.
        all_habitat_idx = np.concatenate(list(habitat_masks.values()))
        hr, hc = np.unravel_index(all_habitat_idx, (rows, cols))
        node_current[hr, hc] = np.nan

        # log10 stretch — Circuitscape convention. Compresses 3+ orders of
        # magnitude into something a linear colour ramp can resolve.
        finite = node_current[np.isfinite(node_current)]
        if finite.size:
            noise_floor = max(float(np.percentile(finite, 1)), 1e-12)
        else:
            noise_floor = 1e-12
        with np.errstate(divide="ignore", invalid="ignore"):
            log_curr = np.log10(np.maximum(node_current, noise_floor))
        log_curr[~np.isfinite(node_current)] = np.nan

        out_arr = np.where(np.isnan(log_curr), -9999.0,
                           log_curr).astype(np.float32)
        out_path = os.path.join(tempfile.gettempdir(),
                                f"wildboar_pinchpoints_{os.getpid()}.tif")
        with rasterio.open(
            out_path, "w",
            driver="GTiff",
            height=rows, width=cols, count=1,
            dtype="float32",
            crs=crs,
            transform=win_tf,
            nodata=-9999,
        ) as dst:
            dst.write(out_arr, 1)
        self.current_raster_path = out_path

    # =================================================================
    # UI thread: layer creation
    # =================================================================
    def finished(self, ok):
        if not ok:
            msg = (f"Failed: {self.exception}" if self.exception
                   else "Cancelled.")
            QMessageBox.critical(None, "Wildboar Connectivity", msg)
            return

        # LCP corridors as one vector layer with attributes per pair.
        if self.lcp_records:
            layer = QgsVectorLayer(
                f"LineString?crs={self.crs_wkt}",
                f"LCP corridors ({len(self.lcp_records)} pairs)",
                "memory")
            prov = layer.dataProvider()
            prov.addAttributes([
                QgsField("from_id", QVariant.Int),
                QgsField("to_id", QVariant.Int),
                QgsField("cost",    QVariant.Double),
            ])
            layer.updateFields()
            for rec in self.lcp_records:
                f = QgsFeature()
                f.setGeometry(QgsGeometry.fromPolylineXY(rec["points"]))
                f.setAttributes([rec["from_id"], rec["to_id"], rec["cost"]])
                prov.addFeature(f)
            layer.updateExtents()
            sym = layer.renderer().symbol()
            sym.setWidth(1.4)
            sym.setColor(QColor("#d62728"))
            QgsProject.instance().addMapLayer(layer)

        # Cumulative current flow / pinchpoints raster.
        if self.current_raster_path:
            rlayer = QgsRasterLayer(
                self.current_raster_path,
                f"Pinchpoints (log10 cum. current, {self.n_pairs} pairs)")
            if rlayer.isValid():
                self._apply_singleband_ramp(rlayer, self.ramp_name)
                QgsProject.instance().addMapLayer(rlayer)

        QgsMessageLog.logMessage(
            f"Done. LCPs={len(self.lcp_records)}, "
            f"pinchpoints={'yes' if self.current_raster_path else 'no'}, "
            f"pairs={self.n_pairs}.", LOG_TAG, Qgis.Success)

    # -----------------------------------------------------------------
    @staticmethod
    def _apply_singleband_ramp(rlayer, ramp_name="Magma"):
        """Pseudocolour renderer with a percentile stretch on the data."""
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
