# -*- coding: utf-8 -*-
"""
WildboarConnectivity v2 - main plugin module.

Workflow
--------
1. Pick a resistance raster and the "Kerneinstaende" habitat polygon
   layer from two dropdowns.
2. Click ONE point on the map. The click defines the centre of a
   circular movement neighbourhood whose radius is the configured
   "wild boar range" (default 7 km; typical home-range radius is 1-4 km,
   dispersal up to ~10 km, so 7 km is a defensible movement window).
3. The plugin selects every habitat polygon intersecting that circle
   and, for every pair of selected habitats, computes:
       - a least-cost path (Dijkstra, scikit-image)
       - a contribution to the cumulative Circuitscape current flow
       - a weighted graph edge (1 / LCP cost) with node centrality

Mathematical conversion of resistance to conductance
----------------------------------------------------
Each pair of 4-connected cells (i, j) becomes one resistor whose
resistance is the average of the two cell resistances (each cell
contributes half the path):

        R_edge(i, j) = 0.5 * (R[i] + R[j])
        g(i, j)      = 1 / R_edge(i, j)

The weighted graph Laplacian L = D - A (with D the diagonal of incident
conductances) is solved against b = +1/n_s on source cells and -1/n_s on
sink cells (the habitat polygons themselves, rasterised onto the grid).
Edge currents are g_ij * (v_i - v_j); node values are half the sum of
|edge currents| on incident edges. Summed across all pairs this gives
the cumulative current map; bright streaks = pinchpoints.
"""

import json
import os
import sys

# ----------------------------------------------------------------------
# Windows / QGIS / NumPy stderr workaround. Must run BEFORE numpy is
# imported anywhere downstream.
# ----------------------------------------------------------------------
if sys.stderr is None:
    class _DummyStderr:
        def write(self, *a, **kw): pass
        def flush(self): pass
    sys.stderr = _DummyStderr()

from qgis.core import (
    Qgis,
    QgsApplication,
    QgsCoordinateTransform,
    QgsFeature,
    QgsFeatureRequest,
    QgsField,
    QgsGeometry,
    QgsMapLayerProxyModel,
    QgsMessageLog,
    QgsPointXY,
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
    QgsWkbTypes,
)
from qgis.PyQt.QtCore import QCoreApplication, QSettings, QTranslator, QVariant
from qgis.PyQt.QtGui import QColor, QIcon
from qgis.PyQt.QtWidgets import QAction, QMessageBox

from .connectivity_dialog import WildboarConnectivityDialog
from .connectivity_task import HabitatConnectivityTask
from .point_tool import PointSelectionTool


LOG_TAG = "wildboar-v2"


# =====================================================================
# Plugin class
# =====================================================================
class WildboarConnectivity:
    """QGIS plugin entry point."""

    def __init__(self, iface):
        self.iface = iface
        self.canvas = iface.mapCanvas()
        self.plugin_dir = os.path.dirname(__file__)

        # Standard i18n scaffold (kept for parity with the original).
        locale = (QSettings().value("locale/userLocale") or "en")[0:2]
        qm = os.path.join(self.plugin_dir, "i18n",
                          f"WildboarConnectivity_{locale}.qm")
        if os.path.exists(qm):
            self.translator = QTranslator()
            self.translator.load(qm)
            QCoreApplication.installTranslator(self.translator)

        self.actions = []
        self.menu = self.tr("&Wildboar Connectivity v2")

        # Runtime state.
        self.dlg = None
        self.point_tool = None
        self.center_point_canvas = None   # QgsPointXY in canvas CRS
        self.active_task = None

    # -----------------------------------------------------------------
    @staticmethod
    def tr(message):
        return QCoreApplication.translate("WildboarConnectivity", message)

    def initGui(self):
        icon_path = os.path.join(self.plugin_dir, "icon.png")
        action = QAction(QIcon(icon_path),
                         self.tr("Habitat connectivity (LCP + pinchpoints + graph)"),
                         self.iface.mainWindow())
        action.triggered.connect(self.run)
        self.iface.addToolBarIcon(action)
        self.iface.addPluginToMenu(self.menu, action)
        self.actions.append(action)

    def unload(self):
        for action in self.actions:
            self.iface.removePluginMenu(self.menu, action)
            self.iface.removeToolBarIcon(action)
        if self.point_tool is not None:
            self.point_tool.clear_marker()

    # -----------------------------------------------------------------
    def run(self):
        """Open the dialog. Wire signals on first invocation only."""
        if self.dlg is None:
            self.dlg = WildboarConnectivityDialog(self.iface.mainWindow())
            self.dlg.cmbRaster.setFilters(QgsMapLayerProxyModel.RasterLayer)
            self.dlg.cmbHabitats.setFilters(QgsMapLayerProxyModel.PolygonLayer)

            self.point_tool = PointSelectionTool(
                self.canvas, role="center", color="#f1c40f")
            self.point_tool.point_selected.connect(self._on_point_selected)

            self.dlg.btnPick.clicked.connect(self._activate_picker)
            self.dlg.btnRun.clicked.connect(self._on_run_clicked)

        self.dlg.show()
        self.dlg.raise_()
        self.dlg.activateWindow()

    # -----------------------------------------------------------------
    # Centre point picking (any point on the map - not tied to a habitat)
    # -----------------------------------------------------------------
    def _activate_picker(self):
        self.dlg.lblStatus.setText(
            "Click anywhere on the map to set the analysis centre...")
        self.point_tool.activate_for_role("center")

    def _on_point_selected(self, point, role):
        self.center_point_canvas = QgsPointXY(point)
        self.dlg.lblRegion.setText(
            f"Centre: {point.x():.1f}, {point.y():.1f}")
        self.dlg.lblStatus.setText(
            "Centre captured. Adjust range and press Run.")

    # -----------------------------------------------------------------
    # Run analysis
    # -----------------------------------------------------------------
    def _on_run_clicked(self):
        raster = self.dlg.cmbRaster.currentLayer()
        habitats = self.dlg.cmbHabitats.currentLayer()

        # ---- Input validation -----------------------------------------
        if raster is None or not isinstance(raster, QgsRasterLayer):
            self._warn("Select a resistance raster.")
            return
        if habitats is None or not isinstance(habitats, QgsVectorLayer):
            self._warn("Select a habitats polygon layer.")
            return
        if habitats.geometryType() != QgsWkbTypes.PolygonGeometry:
            self._warn("Habitat layer must be polygons.")
            return
        if self.center_point_canvas is None:
            self._warn("Pick a centre point on the map.")
            return
        if not (self.dlg.chkLcp.isChecked()
                or self.dlg.chkCircuit.isChecked()
                or self.dlg.chkGraph.isChecked()):
            self._warn("Enable at least one method.")
            return

        # ---- Buffer the click point in the habitats CRS ---------------
        range_m = float(self.dlg.spnRangeKm.value()) * 1000.0
        habitats_crs = habitats.crs()
        if habitats_crs.isGeographic():
            self._warn("Habitats layer is in degrees; reproject to a metric "
                       "CRS (e.g. EPSG:2056) so the range buffer is in meters.")
            return

        canvas_crs = self.canvas.mapSettings().destinationCrs()
        tr_to_habitats = QgsCoordinateTransform(
            canvas_crs, habitats_crs, QgsProject.instance())
        center_in_habitats = tr_to_habitats.transform(self.center_point_canvas)

        center_geom = QgsGeometry.fromPointXY(center_in_habitats)
        # Region = user-defined wild-boar range (drives habitat selection
        # and the visible output extent).
        buffer_geom = center_geom.buffer(range_m, 48)
        # Computation region = region + edge buffer so the circuit solver
        # has resistance cells outside the area of interest, avoiding the
        # boundary artefact where current is artificially squeezed
        # against the window edge. 30 % of the range (min 1.5 km) is
        # typically enough for the radial 1/r decay to settle.
        edge_buffer_m = max(range_m * 0.3, 1500.0)
        compute_buffer_geom = center_geom.buffer(
            range_m + edge_buffer_m, 48)

        # All habitats whose geometry touches the *visible* region.
        # Habitats outside the region (but inside the edge buffer) are
        # NOT sources/sinks - the edge buffer only contributes pure
        # resistance cells to absorb the boundary effect.
        req = (QgsFeatureRequest()
               .setFilterRect(buffer_geom.boundingBox())
               .setFlags(QgsFeatureRequest.ExactIntersect))
        in_range = []
        for f in habitats.getFeatures(req):
            if f.geometry().intersects(buffer_geom):
                in_range.append(f)

        if len(in_range) < 2:
            self._warn(
                f"Only {len(in_range)} habitat(s) found within "
                f"{range_m/1000:.0f} km of the centre. Move the point or "
                f"increase the range.")
            return

        # ---- Transform geometries into the raster CRS ------------------
        raster_crs = raster.crs()
        tr_to_raster = QgsCoordinateTransform(
            habitats_crs, raster_crs, QgsProject.instance())

        # Visible region (in raster CRS): the geojson used to crop the
        # output pinchpoint raster - so its rendered extent matches the
        # user-defined circle exactly.
        buffer_in_raster = QgsGeometry(buffer_geom)
        buffer_in_raster.transform(tr_to_raster)
        buffer_geojson = json.loads(buffer_in_raster.asJson())

        # Computation window (in raster CRS): the bigger buffered circle's
        # bbox. The solver runs on this larger grid; the output is then
        # cropped/masked back to the user region.
        compute_in_raster = QgsGeometry(compute_buffer_geom)
        compute_in_raster.transform(tr_to_raster)
        bbox_in_raster = compute_in_raster.boundingBox()

        habitats_data = []
        for f in in_range:
            g = QgsGeometry(f.geometry())
            g.transform(tr_to_raster)
            geojson = json.loads(g.asJson())
            centroid_pt = g.centroid().asPoint()
            habitats_data.append({
                "id":       int(f.id()),
                "label":    self._describe_feature(habitats, f),
                "geojson":  geojson,
                "area":     float(g.area()),
                "centroid": (float(centroid_pt.x()), float(centroid_pt.y())),
            })
            # Habitats touching but partly outside the circle should
            # still be fully captured by the computational window so
            # their injection cells are not clipped. The output raster
            # is masked back to the circle afterwards.
            bbox_in_raster.combineExtentWith(g.boundingBox())

        if not bbox_in_raster.intersects(raster.extent()):
            self._warn("Habitats do not overlap the raster extent.")
            return

        # ---- Visualize buffer + selected habitats on the map -----------
        if self.dlg.chkShowBuffer.isChecked():
            self._show_region_preview(
                center_in_habitats, buffer_geom, habitats_crs, in_range)

        # ---- Launch background task ------------------------------------
        window_bounds = (bbox_in_raster.xMinimum(),
                         bbox_in_raster.yMinimum(),
                         bbox_in_raster.xMaximum(),
                         bbox_in_raster.yMaximum())

        n_pairs = len(in_range) * (len(in_range) - 1) // 2
        self.dlg.lblStatus.setText(
            f"Running on {len(in_range)} habitats ({n_pairs} pairs)... "
            f"watch the Tasks panel.")

        self.active_task = HabitatConnectivityTask(
            raster_path=raster.source(),
            habitats=habitats_data,
            window_bounds=window_bounds,
            crs_wkt=raster_crs.toWkt(),
            run_lcp=self.dlg.chkLcp.isChecked(),
            run_circuit=self.dlg.chkCircuit.isChecked(),
            run_graph=self.dlg.chkGraph.isChecked(),
            region_geojson=buffer_geojson,
        )
        QgsApplication.taskManager().addTask(self.active_task)

    # -----------------------------------------------------------------
    @staticmethod
    def _describe_feature(layer, feature):
        for name_field in ("name", "Name", "NAME",
                           "bezeichnung", "Bezeichnung",
                           "id", "ID", "fid"):
            i = layer.fields().indexOf(name_field)
            if i >= 0:
                val = feature.attribute(name_field)
                if val not in (None, ""):
                    area_km2 = feature.geometry().area() / 1e6
                    return f"{name_field}={val} (~{area_km2:.2f} km^2)"
        area_km2 = feature.geometry().area() / 1e6
        return f"fid={feature.id()} (~{area_km2:.2f} km^2)"

    # -----------------------------------------------------------------
    def _show_region_preview(self, center_pt, buffer_geom, crs,
                             in_range_features):
        """Add memory layers showing the analysis circle and the
        habitats it selected."""
        authid = crs.authid() or "EPSG:2056"

        # Buffer / range circle.
        buf = QgsVectorLayer(f"Polygon?crs={authid}",
                             "Analysis range (buffer)", "memory")
        prov = buf.dataProvider()
        f = QgsFeature()
        f.setGeometry(buffer_geom)
        prov.addFeature(f)
        buf.updateExtents()
        sym = buf.renderer().symbol()
        sym.setColor(QColor(241, 196, 15, 60))   # translucent yellow
        sym.symbolLayer(0).setStrokeColor(QColor(241, 196, 15))
        QgsProject.instance().addMapLayer(buf)

        # Centre point marker as a permanent layer (cleaner than just
        # the on-canvas vertex marker which disappears after a refresh).
        pt = QgsVectorLayer(f"Point?crs={authid}",
                            "Analysis centre", "memory")
        prov = pt.dataProvider()
        ff = QgsFeature()
        ff.setGeometry(QgsGeometry.fromPointXY(center_pt))
        prov.addFeature(ff)
        pt.updateExtents()
        sym = pt.renderer().symbol()
        sym.setColor(QColor(241, 196, 15))
        sym.setSize(4.0)
        QgsProject.instance().addMapLayer(pt)

        # Selected habitats highlighted in green.
        sel = QgsVectorLayer(
            f"Polygon?crs={authid}",
            f"Selected habitats ({len(in_range_features)})", "memory")
        prov = sel.dataProvider()
        prov.addAttributes([QgsField("id", QVariant.Int)])
        sel.updateFields()
        for fe in in_range_features:
            ff = QgsFeature()
            ff.setGeometry(fe.geometry())
            ff.setAttributes([fe.id()])
            prov.addFeature(ff)
        sel.updateExtents()
        sym = sel.renderer().symbol()
        sym.setColor(QColor(46, 204, 113, 100))  # translucent green
        QgsProject.instance().addMapLayer(sel)

    # -----------------------------------------------------------------
    def _warn(self, msg):
        self.dlg.lblStatus.setText(msg)
        QgsMessageLog.logMessage(msg, LOG_TAG, Qgis.Warning)
        QMessageBox.warning(self.dlg, "Wildboar Connectivity", msg)
