# -*- coding: utf-8 -*-
"""
WildboarConnectivity v2 - main plugin module.

Workflow
--------
1. User picks a resistance raster and the "Kerneinstaende" habitat polygon
   layer from two dropdowns.
2. User clicks one point on the canvas. The plugin selects the habitat
   polygon that contains (or is nearest to) that click.
3. The plugin buffers the selected habitat by a configurable "wild boar
   range" (default 5 km) and collects every habitat polygon that
   intersects the buffer. This restricts the analysis to a spatial scale
   that is ecologically defensible for wild boar (typical daily movement
   1-5 km, mean dispersal ~10 km).
4. For all pairs of habitats inside the buffer the plugin computes:
       - least-cost paths (Dijkstra, scikit-image)
       - cumulative Circuitscape-style current flow (pinchpoints)

Mathematical conversion of resistance to conductance
----------------------------------------------------
Each pair of 4-connected cells (i, j) becomes one resistor whose
resistance is the average of the two cell resistances (each cell
contributes half the path):

        R_edge(i, j) = 0.5 * (R[i] + R[j])
        g(i, j)      = 1 / R_edge(i, j)

The weighted graph Laplacian L = D - A (with D the diagonal of incident
conductances) is solved against b = +1/n_s on source cells and -1/n_s on
sink cells. Edge currents are g_ij * (v_i - v_j); node values are half
the sum of |edge currents| on incident edges. The result is summed
across all source/sink pairs to give the *cumulative* current map.
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
        self.target_feature_id = None    # fid of the clicked habitat
        self.target_layer_id = None      # which layer it was clicked on
        self.active_task = None

    # -----------------------------------------------------------------
    @staticmethod
    def tr(message):
        return QCoreApplication.translate("WildboarConnectivity", message)

    def initGui(self):
        icon_path = os.path.join(self.plugin_dir, "icon.png")
        action = QAction(QIcon(icon_path),
                         self.tr("Habitat connectivity (LCP + pinchpoints)"),
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
                self.canvas, role="region", color="#f1c40f")
            self.point_tool.point_selected.connect(self._on_point_selected)

            self.dlg.btnPick.clicked.connect(self._activate_picker)
            self.dlg.btnRun.clicked.connect(self._on_run_clicked)

        self.dlg.show()
        self.dlg.raise_()
        self.dlg.activateWindow()

    # -----------------------------------------------------------------
    # Habitat picking
    # -----------------------------------------------------------------
    def _activate_picker(self):
        if self.dlg.cmbHabitats.currentLayer() is None:
            self._warn("Pick the habitats polygon layer first.")
            return
        self.dlg.lblStatus.setText("Click a habitat on the map...")
        self.point_tool.activate_for_role("region")

    def _on_point_selected(self, point, role):
        habitats = self.dlg.cmbHabitats.currentLayer()
        if habitats is None:
            self._warn("No habitats layer selected.")
            return

        # Transform canvas coords -> habitats CRS for the spatial query.
        canvas_crs = self.canvas.mapSettings().destinationCrs()
        tr = QgsCoordinateTransform(canvas_crs, habitats.crs(),
                                    QgsProject.instance())
        pt = tr.transform(point)
        target = self._find_target_feature(habitats, pt)

        if target is None:
            self._warn("No habitat found near the click.")
            return

        self.target_feature_id = target.id()
        self.target_layer_id = habitats.id()
        self.dlg.lblRegion.setText(
            f"Selected: {self._describe_feature(habitats, target)}")
        self.dlg.lblStatus.setText(
            "Habitat selected. Adjust range and press Run.")

    @staticmethod
    def _find_target_feature(layer, point_xy):
        """Return the polygon containing the point, else the nearest."""
        pt_geom = QgsGeometry.fromPointXY(point_xy)
        nearest = None
        nearest_d = float("inf")
        # Pre-filter by bbox + small expansion to speed things up on big layers.
        for f in layer.getFeatures():
            g = f.geometry()
            if g is None or g.isEmpty():
                continue
            if g.contains(pt_geom):
                return f
            d = g.distance(pt_geom)
            if d < nearest_d:
                nearest_d = d
                nearest = f
        return nearest

    @staticmethod
    def _describe_feature(layer, feature):
        """Best-effort human-readable name for a habitat feature."""
        for name_field in ("name", "Name", "NAME",
                           "bezeichnung", "Bezeichnung",
                           "id", "ID", "fid"):
            i = layer.fields().indexOf(name_field)
            if i >= 0:
                val = feature.attribute(name_field)
                if val not in (None, ""):
                    area_km2 = feature.geometry().area() / 1e6
                    return f"{name_field}={val}  (~{area_km2:.2f} km^2)"
        area_km2 = feature.geometry().area() / 1e6
        return f"fid={feature.id()}  (~{area_km2:.2f} km^2)"

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
        if (habitats.geometryType() != QgsWkbTypes.PolygonGeometry):
            self._warn("Habitat layer must be polygons.")
            return
        if self.target_feature_id is None:
            self._warn("Pick a starting habitat on the map.")
            return
        if not (self.dlg.chkLcp.isChecked()
                or self.dlg.chkCircuit.isChecked()):
            self._warn("Enable at least one method.")
            return

        target = habitats.getFeature(self.target_feature_id)
        if not target.isValid():
            self._warn("Selected habitat no longer exists in the layer.")
            return

        # ---- Build buffer in habitat CRS -------------------------------
        range_m = float(self.dlg.spnRangeKm.value()) * 1000.0
        habitats_crs = habitats.crs()
        if habitats_crs.isGeographic():
            self._warn("Habitats layer is in degrees; reproject to a metric "
                       "CRS (e.g. EPSG:2056) so the range buffer is in meters.")
            return

        buffer_geom = target.geometry().buffer(range_m, 24)

        # All habitat features whose geometry touches the buffer.
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
                f"{range_m/1000:.0f} km. Increase the range.")
            return

        # ---- Transform geometries into the raster CRS ------------------
        raster_crs = raster.crs()
        tr = QgsCoordinateTransform(habitats_crs, raster_crs,
                                    QgsProject.instance())

        habitats_data = []
        bbox_in_raster = None
        for f in in_range:
            g = QgsGeometry(f.geometry())
            g.transform(tr)
            geojson = json.loads(g.asJson())
            habitats_data.append({
                "id": int(f.id()),
                "label": self._describe_feature(habitats, f),
                "geojson": geojson,
                "area": float(g.area()),
            })
            bb = g.boundingBox()
            if bbox_in_raster is None:
                bbox_in_raster = bb
            else:
                bbox_in_raster.combineExtentWith(bb)

        # Pad the analysis window so corridors can bend outside the
        # tight habitat envelope. 500 m is generous for wild boar.
        bbox_in_raster.grow(500.0)

        # Confirm the window overlaps the raster.
        if not bbox_in_raster.intersects(raster.extent()):
            self._warn("Habitats do not overlap the raster extent.")
            return

        # ---- Visualize buffer + selected habitats on the map -----------
        if self.dlg.chkShowBuffer.isChecked():
            self._show_region_preview(buffer_geom, habitats_crs, in_range)

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
            target_id=int(self.target_feature_id),
            window_bounds=window_bounds,
            crs_wkt=raster_crs.toWkt(),
            run_lcp=self.dlg.chkLcp.isChecked(),
            run_circuit=self.dlg.chkCircuit.isChecked(),
        )
        QgsApplication.taskManager().addTask(self.active_task)

    # -----------------------------------------------------------------
    def _show_region_preview(self, buffer_geom, crs, in_range_features):
        """Add two memory layers showing what region is being analyzed."""
        authid = crs.authid() or "EPSG:2056"

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
