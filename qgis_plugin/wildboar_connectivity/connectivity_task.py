# -*- coding: utf-8 -*-
"""
Background QgsTask that computes - for every pair of habitats inside the
user-defined wild boar movement neighbourhood:

    * Pairwise least-cost paths (Dijkstra, scikit-image).
    * Cumulative Circuitscape-style current-flow / pinchpoint raster.
    * A weighted habitat NETWORK GRAPH: habitat centroids as nodes,
      straight edges weighted by 1 / LCP cost, with degree and
      betweenness-centrality scores per node (networkx).

Math summary (full derivation in main_plugin.py docstring):
    g_ij = 1 / (0.5 * (R_i + R_j))
    L    = D - A  (weighted graph Laplacian, D = sum of incident g)
    L v = b   with  b[A] = +1/|A|, b[B] = -1/|B| on habitat polygons A, B
    I_e  = g_ij * (v_i - v_j)
    Cum_e = sum over habitat pairs of |I_e|.
    Pixel value = 0.5 * sum_{j ~ i} Cum(i,j).

Performance: the Laplacian is factored ONCE (scipy splu) with a single
corner cell grounded; every pair is then a fast back-substitution.
"""

import os
import tempfile
import itertools

import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.windows import from_bounds, Window, transform as window_transform
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
    """Pairwise LCPs, cumulative current flow, and a network graph
    across every habitat inside the analysis window."""

    def __init__(self,
                 raster_path: str,
                 habitats: list,         # [{id, label, geojson, area, centroid}]
                 window_bounds: tuple,   # (min_x, min_y, max_x, max_y) raster CRS
                 crs_wkt: str,           # raster CRS for output layers
                 run_lcp: bool,
                 run_circuit: bool,
                 run_graph: bool,
                 region_geojson: dict = None,   # buffer polygon in raster CRS
                 ramp_name: str = "Magma"):
        super().__init__("Wildboar habitat connectivity", QgsTask.CanCancel)
        self.raster_path = raster_path
        self.habitats = habitats
        self.window_bounds = window_bounds
        self.crs_wkt = crs_wkt
        self.run_lcp = run_lcp
        self.run_circuit = run_circuit
        self.run_graph = run_graph
        self.region_geojson = region_geojson
        self.ramp_name = ramp_name

        # Outputs filled by run(), consumed by finished().
        self.lcp_records = []           # [{from_id, to_id, cost, points}]
        self.current_raster_path = None
        self.graph_nodes = []           # [{id, centroid_xy, degree, betweenness, area}]
        self.graph_edges = []           # [{from_id, to_id, cost, weight, from_xy, to_xy}]
        self.n_pairs = 0
        self.n_habitats = 0
        self.exception = None

    # =================================================================
    # Worker thread
    # =================================================================
    def run(self):
        try:
            arr, win_tf, crs, rows, cols, _nodata_mask = self._load_window()

            habitat_masks = self._rasterize_habitats(arr.shape, win_tf)
            self.n_habitats = len(habitat_masks)
            if self.n_habitats < 2:
                raise RuntimeError(
                    f"Only {self.n_habitats} habitats rasterise onto the "
                    f"raster window. Need at least 2.")

            pairs = list(itertools.combinations(habitat_masks.keys(), 2))
            self.n_pairs = len(pairs)
            QgsMessageLog.logMessage(
                f"Analysing {self.n_habitats} habitats, {self.n_pairs} pairs "
                f"on a {rows}x{cols} window.", LOG_TAG, Qgis.Info)

            # ---- LCPs (cheap; also used by the graph step) -----------
            need_lcp = self.run_lcp or self.run_graph
            lcp_costs = {}
            if need_lcp:
                lcp_costs = self._compute_lcps(arr, win_tf, habitat_masks, pairs)
                if not self.run_lcp:
                    # User wanted graph only - don't add LCP layer, but
                    # keep costs for the graph step.
                    self.lcp_records = []
            self.setProgress(35)
            if self.isCanceled():
                return False

            # ---- Cumulative current flow ------------------------------
            if self.run_circuit:
                self._compute_cumulative_currents(
                    arr, win_tf, crs, habitat_masks, pairs)
            self.setProgress(85)
            if self.isCanceled():
                return False

            # ---- Network graph (needs lcp_costs) ----------------------
            if self.run_graph:
                self._build_graph(habitat_masks, pairs, lcp_costs)
            self.setProgress(100)
            return True

        except Exception as exc:
            self.exception = exc
            QgsMessageLog.logMessage(f"Task failed: {exc}",
                                     LOG_TAG, Qgis.Critical)
            return False

    # -----------------------------------------------------------------
    def _load_window(self):
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
        arr = np.clip(arr, 1e-3, 1e6)
        return arr, win_tf, crs, arr.shape[0], arr.shape[1], nodata_mask

    def _rasterize_habitats(self, shape, win_tf):
        out = {}
        for h in self.habitats:
            mask = rasterize(
                [(h["geojson"], 1)],
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

    def _habitat_anchor_rc(self, flat_ids, rows, cols):
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
        """Run Dijkstra for every pair. Returns {(a,b): cost} regardless
        of whether the LCP layer itself will be added."""
        from skimage.graph import route_through_array
        rows, cols = arr.shape

        anchors = {hid: self._habitat_anchor_rc(idx, rows, cols)
                   for hid, idx in habitat_masks.items()}

        lcp_costs = {}
        for k, (a, b) in enumerate(pairs):
            if self.isCanceled():
                return lcp_costs
            try:
                indices, cost = route_through_array(
                    arr, anchors[a], anchors[b],
                    fully_connected=True, geometric=True)
                pts = [QgsPointXY(*rio_xy(win_tf, r, c)) for r, c in indices]
                lcp_costs[(a, b)] = float(cost)
                self.lcp_records.append({
                    "from_id": int(a), "to_id": int(b),
                    "cost": float(cost), "points": pts,
                })
            except Exception as exc:
                QgsMessageLog.logMessage(
                    f"LCP {a}<->{b} failed: {exc}", LOG_TAG, Qgis.Warning)
            self.setProgress(5 + 25 * (k + 1) / max(1, self.n_pairs))
        return lcp_costs

    # -----------------------------------------------------------------
    def _compute_cumulative_currents(self, arr, win_tf, crs,
                                     habitat_masks, pairs):
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import splu

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

        # Ground one cell and factor ONCE.
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
            bvec[a_ids] += 1.0 / a_ids.size
            bvec[b_ids] -= 1.0 / b_ids.size
            b_red = bvec[keep]

            v_red = lu.solve(b_red)
            v = np.zeros(n)
            v[keep] = v_red

            edge_I = g_edge * (v[e_i] - v[e_j])
            cum_edge += np.abs(edge_I)
            self.setProgress(55 + 25 * (k + 1) / max(1, self.n_pairs))

        node_current = np.zeros(n)
        np.add.at(node_current, e_i, cum_edge)
        np.add.at(node_current, e_j, cum_edge)
        node_current *= 0.5
        node_current = node_current.reshape(rows, cols)

        # Mask habitat cells (core habitat, not corridor).
        all_habitat_idx = np.concatenate(list(habitat_masks.values()))
        hr, hc = np.unravel_index(all_habitat_idx, (rows, cols))
        node_current[hr, hc] = np.nan

        # ---- Clip to the user-defined circular region ----------------
        # Computation ran on the bigger habitat-bbox window (so partial
        # habitats are not clipped); the OUTPUT raster is masked to the
        # circle so its visible extent matches the analysis range
        # exactly.
        if self.region_geojson is not None:
            inside_region = rasterize(
                [(self.region_geojson, 1)],
                out_shape=node_current.shape,
                transform=win_tf,
                fill=0,
                dtype="uint8",
                all_touched=False,
            ).astype(bool)
            node_current[~inside_region] = np.nan
        else:
            inside_region = np.ones(node_current.shape, dtype=bool)

        # log10 stretch (Circuitscape convention).
        finite = node_current[np.isfinite(node_current)]
        if finite.size:
            noise_floor = max(float(np.percentile(finite, 1)), 1e-12)
        else:
            noise_floor = 1e-12
        with np.errstate(divide="ignore", invalid="ignore"):
            log_curr = np.log10(np.maximum(node_current, noise_floor))
        log_curr[~np.isfinite(node_current)] = np.nan

        # ---- Crop the GeoTIFF tightly around the circle --------------
        # Find the bounding box of cells inside the region; the output
        # raster's extent ends up as the minimum rectangle that still
        # contains the full circle (no surrounding empty padding).
        valid_rows = np.where(inside_region.any(axis=1))[0]
        valid_cols = np.where(inside_region.any(axis=0))[0]
        if valid_rows.size and valid_cols.size:
            r0, r1 = int(valid_rows[0]), int(valid_rows[-1]) + 1
            c0, c1 = int(valid_cols[0]), int(valid_cols[-1]) + 1
        else:
            r0, r1, c0, c1 = 0, rows, 0, cols
        log_curr = log_curr[r0:r1, c0:c1]
        out_tf = window_transform(Window(c0, r0, c1 - c0, r1 - r0), win_tf)
        out_h, out_w = log_curr.shape

        out_arr = np.where(np.isnan(log_curr), -9999.0,
                           log_curr).astype(np.float32)
        out_path = os.path.join(tempfile.gettempdir(),
                                f"wildboar_pinchpoints_{os.getpid()}.tif")
        with rasterio.open(
            out_path, "w",
            driver="GTiff",
            height=out_h, width=out_w, count=1,
            dtype="float32",
            crs=crs,
            transform=out_tf,
            nodata=-9999,
        ) as dst:
            dst.write(out_arr, 1)
        self.current_raster_path = out_path

    # -----------------------------------------------------------------
    def _build_graph(self, habitat_masks, pairs, lcp_costs):
        """Build a NetworkX graph with habitats as nodes and 1/LCP cost
        as edge weight. Compute degree and betweenness centrality.

        Why these metrics?
            * degree           = how many other habitats this one connects
                                 to (with finite cost). Robustness indicator.
            * betweenness      = fraction of all habitat-to-habitat
                                 shortest paths that pass through this
                                 habitat. High = critical stepping stone.

        Outputs are stored as plain dicts; the UI thread builds the
        actual QgsVectorLayers in finished().
        """
        try:
            import networkx as nx
        except ImportError:
            QgsMessageLog.logMessage(
                "networkx not installed; skipping graph step.",
                LOG_TAG, Qgis.Warning)
            return

        # Look up centroids from the habitats payload we were given.
        centroid_by_id = {h["id"]: h["centroid"] for h in self.habitats
                          if h["id"] in habitat_masks}
        area_by_id = {h["id"]: h["area"] for h in self.habitats
                      if h["id"] in habitat_masks}

        G = nx.Graph()
        for hid in habitat_masks.keys():
            G.add_node(hid)

        for (a, b) in pairs:
            cost = lcp_costs.get((a, b))
            if cost is None or cost <= 0 or not np.isfinite(cost):
                continue
            G.add_edge(a, b,
                       cost=cost,
                       weight=1.0 / cost)  # higher = stronger connection

        # Shortest-path distance for betweenness uses raw LCP cost.
        try:
            betw = nx.betweenness_centrality(G, weight="cost", normalized=True)
        except Exception:
            betw = {n: 0.0 for n in G.nodes}
        degree = dict(G.degree())

        for hid in G.nodes:
            self.graph_nodes.append({
                "id":          int(hid),
                "centroid_xy": centroid_by_id.get(hid, (0.0, 0.0)),
                "degree":      int(degree.get(hid, 0)),
                "betweenness": float(betw.get(hid, 0.0)),
                "area":        float(area_by_id.get(hid, 0.0)),
            })

        for (a, b, d) in G.edges(data=True):
            self.graph_edges.append({
                "from_id": int(a),
                "to_id":   int(b),
                "cost":    float(d["cost"]),
                "weight":  float(d["weight"]),
                "from_xy": centroid_by_id.get(a, (0.0, 0.0)),
                "to_xy":   centroid_by_id.get(b, (0.0, 0.0)),
            })

    # =================================================================
    # UI thread: layer creation
    # =================================================================
    def finished(self, ok):
        if not ok:
            msg = (f"Failed: {self.exception}" if self.exception
                   else "Cancelled.")
            QMessageBox.critical(None, "Wildboar Connectivity", msg)
            return

        # --- LCP corridor layer ----------------------------------------
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

        # --- Cumulative current / pinchpoints raster -------------------
        if self.current_raster_path:
            rlayer = QgsRasterLayer(
                self.current_raster_path,
                f"Pinchpoints (log10 cum. current, {self.n_pairs} pairs)")
            if rlayer.isValid():
                self._apply_singleband_ramp(rlayer, self.ramp_name)
                QgsProject.instance().addMapLayer(rlayer)

        # --- Graph: edges & nodes --------------------------------------
        if self.graph_edges:
            self._add_graph_edges_layer()
        if self.graph_nodes:
            self._add_graph_nodes_layer()

        QgsMessageLog.logMessage(
            f"Done. habitats={self.n_habitats}, pairs={self.n_pairs}, "
            f"LCPs={len(self.lcp_records)}, "
            f"pinchpoints={'yes' if self.current_raster_path else 'no'}, "
            f"graph={'yes' if self.graph_edges else 'no'}.",
            LOG_TAG, Qgis.Success)

    # -----------------------------------------------------------------
    def _add_graph_edges_layer(self):
        layer = QgsVectorLayer(
            f"LineString?crs={self.crs_wkt}",
            f"Habitat graph - edges ({len(self.graph_edges)})",
            "memory")
        prov = layer.dataProvider()
        prov.addAttributes([
            QgsField("from_id", QVariant.Int),
            QgsField("to_id",   QVariant.Int),
            QgsField("cost",    QVariant.Double),
            QgsField("weight",  QVariant.Double),
        ])
        layer.updateFields()
        for e in self.graph_edges:
            f = QgsFeature()
            f.setGeometry(QgsGeometry.fromPolylineXY([
                QgsPointXY(*e["from_xy"]),
                QgsPointXY(*e["to_xy"]),
            ]))
            f.setAttributes([e["from_id"], e["to_id"], e["cost"], e["weight"]])
            prov.addFeature(f)
        layer.updateExtents()
        sym = layer.renderer().symbol()
        sym.setColor(QColor(52, 152, 219, 180))   # translucent blue
        sym.setWidth(0.6)
        QgsProject.instance().addMapLayer(layer)

    def _add_graph_nodes_layer(self):
        layer = QgsVectorLayer(
            f"Point?crs={self.crs_wkt}",
            f"Habitat graph - nodes ({len(self.graph_nodes)})",
            "memory")
        prov = layer.dataProvider()
        prov.addAttributes([
            QgsField("id",          QVariant.Int),
            QgsField("degree",      QVariant.Int),
            QgsField("betweenness", QVariant.Double),
            QgsField("area_km2",    QVariant.Double),
        ])
        layer.updateFields()
        for n in self.graph_nodes:
            f = QgsFeature()
            x, y = n["centroid_xy"]
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
            f.setAttributes([
                n["id"], n["degree"], n["betweenness"], n["area"] / 1e6,
            ])
            prov.addFeature(f)
        layer.updateExtents()
        sym = layer.renderer().symbol()
        sym.setColor(QColor("#2c3e50"))
        sym.setSize(3.5)
        QgsProject.instance().addMapLayer(layer)

    # -----------------------------------------------------------------
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
