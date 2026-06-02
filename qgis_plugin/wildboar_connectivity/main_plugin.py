# -*- coding: utf-8 -*-
"""
WildboarConnectivity (ASF) - main plugin module.

Habitat-free decision-support tool for African Swine Fever:

    1.  Load a resistance raster.
    2.  Click ONE point on the map. That's the outbreak origin.
    3.  Draw fences and/or place wildlife overpasses (a working copy of
        the resistance grid is modified, never the source GeoTIFF).
    4.  Press Run. The plugin produces:
            - Pinchpoints / corridors raster (Circuit theory: origin
              disc as source, AOI boundary as sink) - the ASF spread
              bottlenecks; bright = best fence target.
            - Continuous infection-risk raster (green / yellow / red).
            - Optional: cost-from-origin raster, random walk traces.
    5.  Tweak fences / overpasses, hit Run again, compare.
"""

import json
import os
import sys

# Windows / QGIS / NumPy stderr workaround. Must run BEFORE numpy.
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
    QgsGeometry,
    QgsMapLayerProxyModel,
    QgsMessageLog,
    QgsPointXY,
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
    QgsWkbTypes,
)
from qgis.PyQt.QtCore import QCoreApplication, QSettings, QTranslator
from qgis.PyQt.QtGui import QColor, QIcon
from qgis.PyQt.QtWidgets import QAction, QMessageBox

from .connectivity_dialog import WildboarConnectivityDialog
from .analysis.connectivity_task import AsfConnectivityTask
from .tools.fence_tool import FenceDrawingTool
from .utils.layer_utils import (
    ZOrder,
    add_wildboar_layer,
    ensure_swissimage,
    remove_all_wildboar_layers,
)
from .tools.point_tool import PointSelectionTool


LOG_TAG = "wildboar-asf"

# Generous initial read window around the origin (metres). The task
# computes the actual AOI from landscape cost, so this number only sets
# the maximum extent the analysis could ever look at.
INITIAL_BUFFER_M = 10000.0
EDGE_BUFFER_M    = 1500.0


class WildboarConnectivity:
    """QGIS plugin entry point."""

    def __init__(self, iface):
        self.iface = iface
        self.canvas = iface.mapCanvas()
        self.plugin_dir = os.path.dirname(__file__)

        # i18n scaffold.
        locale = (QSettings().value("locale/userLocale") or "en")[0:2]
        qm = os.path.join(self.plugin_dir, "i18n",
                          f"WildboarConnectivity_{locale}.qm")
        if os.path.exists(qm):
            self.translator = QTranslator()
            self.translator.load(qm)
            QCoreApplication.installTranslator(self.translator)

        self.actions = []
        self.menu = self.tr("&Wildboar Connectivity (ASF)")

        # State
        self.dlg = None
        self.origin_tool = None       # PointSelectionTool, role="origin"
        self.overpass_tool = None     # PointSelectionTool, role="overpass"
        self.fence_tool = None        # FenceDrawingTool
        self.active_task = None

        self.origin_canvas_point = None   # QgsPointXY in canvas CRS
        # Fences / overpasses are stored in the canvas CRS and transformed
        # to the raster CRS just before launching the task.
        self.fences = []                  # list[QgsGeometry] polyline
        self.overpasses = []              # list[QgsPointXY]

    # -----------------------------------------------------------------
    @staticmethod
    def tr(msg):
        return QCoreApplication.translate("WildboarConnectivity", msg)

    def initGui(self):
        icon_path = os.path.join(self.plugin_dir, "icon.png")
        action = QAction(
            QIcon(icon_path),
            self.tr("ASF connectivity (no habitats needed)"),
            self.iface.mainWindow())
        action.triggered.connect(self.run)
        self.iface.addToolBarIcon(action)
        self.iface.addPluginToMenu(self.menu, action)
        self.actions.append(action)

    def unload(self):
        for action in self.actions:
            self.iface.removePluginMenu(self.menu, action)
            self.iface.removeToolBarIcon(action)
        for tool in (self.origin_tool, self.overpass_tool, self.fence_tool):
            if tool is not None:
                try:
                    tool.clear_marker() if hasattr(tool, "clear_marker") \
                        else tool.reset()
                except Exception:
                    pass
        try:
            remove_all_wildboar_layers(keep_basemap=False)
        except Exception:
            pass

    # -----------------------------------------------------------------
    def run(self):
        if self.dlg is None:
            self.dlg = WildboarConnectivityDialog(self.iface.mainWindow())
            self.dlg.cmbRaster.setFilters(QgsMapLayerProxyModel.RasterLayer)
            # Habitat polygons are optional (LCP destinations override).
            self.dlg.cmbHabitats.setFilters(QgsMapLayerProxyModel.PolygonLayer)
            self.dlg.cmbHabitats.setAllowEmptyLayer(True)

            self.origin_tool = PointSelectionTool(
                self.canvas, role="origin", color="#e74c3c")
            self.origin_tool.point_selected.connect(self._on_origin_selected)

            self.overpass_tool = PointSelectionTool(
                self.canvas, role="overpass", color="#3498db")
            self.overpass_tool.point_selected.connect(
                self._on_overpass_selected)

            self.fence_tool = FenceDrawingTool(self.canvas, color="#e74c3c")
            self.fence_tool.fence_completed.connect(self._on_fence_drawn)

            self.dlg.btnPickOrigin.clicked.connect(
                lambda: self._activate_tool("origin"))
            self.dlg.btnDrawFence.clicked.connect(
                lambda: self._activate_tool("fence"))
            self.dlg.btnPlaceOverpass.clicked.connect(
                lambda: self._activate_tool("overpass"))
            self.dlg.btnResetMods.clicked.connect(self._on_reset_mods)
            self.dlg.btnRun.clicked.connect(self._on_run_clicked)

        self.dlg.show()
        self.dlg.raise_()
        self.dlg.activateWindow()

    # -----------------------------------------------------------------
    # Tool activation
    # -----------------------------------------------------------------
    def _activate_tool(self, which):
        if which == "origin":
            self.dlg.lblStatus.setText("Click on the map - that is the outbreak origin.")
            self.origin_tool.activate_for_role("origin")
        elif which == "fence":
            self.dlg.lblStatus.setText(
                "Draw fence: left-click to add vertices, "
                "right-click or double-click to finish.")
            self.fence_tool.activate_for_drawing()
        elif which == "overpass":
            self.dlg.lblStatus.setText(
                "Click on the map to place a wildlife overpass.")
            self.overpass_tool.activate_for_role("overpass")

    def _on_origin_selected(self, point_canvas, role):
        self.origin_canvas_point = QgsPointXY(point_canvas)
        self.dlg.lblOrigin.setText(
            f"Origin: {point_canvas.x():.1f}, {point_canvas.y():.1f}")
        self.dlg.lblStatus.setText("Origin captured. Press Run.")

    def _on_fence_drawn(self, geom_canvas):
        self.fences.append(QgsGeometry(geom_canvas))
        self._refresh_mods_label()
        self._refresh_mods_preview_layer()
        self.dlg.lblStatus.setText(
            "Fence added. Press Run to recompute, or draw another.")

    def _on_overpass_selected(self, point_canvas, role):
        self.overpasses.append(QgsPointXY(point_canvas))
        self._refresh_mods_label()
        self._refresh_mods_preview_layer()
        self.dlg.lblStatus.setText(
            "Overpass added. Press Run to recompute, or place another.")

    def _on_reset_mods(self):
        self.fences = []
        self.overpasses = []
        self._refresh_mods_label()
        self._refresh_mods_preview_layer()
        self.dlg.lblStatus.setText("Modifications cleared.")

    def _refresh_mods_label(self):
        self.dlg.lblMods.setText(
            f"Active: {len(self.fences)} fences, "
            f"{len(self.overpasses)} overpasses.")

    # -----------------------------------------------------------------
    def _refresh_mods_preview_layer(self):
        """Live preview of fences and overpasses on the canvas."""
        raster = self.dlg.cmbRaster.currentLayer()
        if raster is None:
            return
        crs_authid = self.canvas.mapSettings().destinationCrs().authid() \
                     or "EPSG:2056"

        # Wipe previous previews.
        proj = QgsProject.instance()
        for lid, lyr in list(proj.mapLayers().items()):
            if lyr.name() in ("Fences (active)", "Overpasses (active)"):
                proj.removeMapLayer(lid)

        if self.fences:
            lay = QgsVectorLayer(
                f"LineString?crs={crs_authid}",
                "Fences (active)", "memory")
            prov = lay.dataProvider()
            for g in self.fences:
                f = QgsFeature()
                f.setGeometry(g)
                prov.addFeature(f)
            lay.updateExtents()
            sym = lay.renderer().symbol()
            sym.setColor(QColor("#e74c3c"))
            sym.setWidth(2.0)
            add_wildboar_layer(lay, ZOrder.FENCES)

        if self.overpasses:
            lay = QgsVectorLayer(
                f"Point?crs={crs_authid}",
                "Overpasses (active)", "memory")
            prov = lay.dataProvider()
            for p in self.overpasses:
                f = QgsFeature()
                f.setGeometry(QgsGeometry.fromPointXY(p))
                prov.addFeature(f)
            lay.updateExtents()
            sym = lay.renderer().symbol()
            sym.setColor(QColor("#3498db"))
            sym.setSize(3.5)
            add_wildboar_layer(lay, ZOrder.FENCES)

    # -----------------------------------------------------------------
    # Run
    # -----------------------------------------------------------------
    def _on_run_clicked(self):
        raster = self.dlg.cmbRaster.currentLayer()

        if raster is None or not isinstance(raster, QgsRasterLayer):
            self._warn("Select a resistance raster.")
            return
        if self.origin_canvas_point is None:
            self._warn("Pick an outbreak point on the map.")
            return
        if not (self.dlg.chkLcp.isChecked()
                or self.dlg.chkCircuit.isChecked()
                or self.dlg.chkRisk.isChecked()
                or self.dlg.chkRandomWalk.isChecked()
                or self.dlg.chkIBMM.isChecked()):
            self._warn("Enable at least one output.")
            return

        canvas_crs = self.canvas.mapSettings().destinationCrs()
        raster_crs = raster.crs()
        if raster_crs.isGeographic():
            self._warn("Resistance raster is in degrees - reproject it to "
                       "a metric CRS (e.g. EPSG:2056).")
            return

        # ---- Transform the origin into raster CRS --------------------
        tr_to_raster = QgsCoordinateTransform(
            canvas_crs, raster_crs, QgsProject.instance())
        origin_in_raster = tr_to_raster.transform(self.origin_canvas_point)
        origin_xy = (float(origin_in_raster.x()),
                     float(origin_in_raster.y()))

        # ---- Initial read window: 10 km around origin ----------------
        ext = raster.extent()
        min_x = max(ext.xMinimum(),
                    origin_xy[0] - (INITIAL_BUFFER_M + EDGE_BUFFER_M))
        max_x = min(ext.xMaximum(),
                    origin_xy[0] + (INITIAL_BUFFER_M + EDGE_BUFFER_M))
        min_y = max(ext.yMinimum(),
                    origin_xy[1] - (INITIAL_BUFFER_M + EDGE_BUFFER_M))
        max_y = min(ext.yMaximum(),
                    origin_xy[1] + (INITIAL_BUFFER_M + EDGE_BUFFER_M))
        if max_x <= min_x or max_y <= min_y:
            self._warn("Outbreak point is outside the resistance raster.")
            return

        # ---- Transform fences / overpasses into raster CRS -----------
        fences_geojson = []
        for fg in self.fences:
            gg = QgsGeometry(fg)
            gg.transform(tr_to_raster)
            fences_geojson.append(json.loads(gg.asJson()))

        overpasses_xy = []
        for p in self.overpasses:
            pp = tr_to_raster.transform(p)
            overpasses_xy.append((float(pp.x()), float(pp.y())))

        # ---- Optional habitat layer for LCP destinations -------------
        # When the user picks a polygon layer in the dialog, the LCP
        # step targets those polygons instead of auto-detecting
        # low-resistance clusters from the raster.
        habitats_geojson = self._collect_habitat_polygons(tr_to_raster)

        # ---- Clean previous outputs, keep Swissimage basemap --------
        remove_all_wildboar_layers(keep_basemap=True)
        ensure_swissimage(opacity=0.6)
        self._show_origin_preview(origin_in_raster, raster_crs)
        self._refresh_mods_preview_layer()

        # ---- Launch background task ----------------------------------
        options = {
            "lcp":                 self.dlg.chkLcp.isChecked(),
            "n_lcp_targets":       100,
            "circuit":             self.dlg.chkCircuit.isChecked(),
            "use_circuitscape_jl": True,   # always, never user-controlled
            "risk":                self.dlg.chkRisk.isChecked(),
            "random_walk":         self.dlg.chkRandomWalk.isChecked(),
            "n_walks":             int(self.dlg.spnNWalks.value()),
            "walk_kappa":          2.0,    # moderate directional persistence
            "ibmm":                self.dlg.chkIBMM.isChecked(),
            "n_agents":            int(self.dlg.spnNAgents.value()),   # default 2000
            "habitats_geojson":    habitats_geojson,
        }
        self.dlg.lblStatus.setText(
            f"Running... {len(self.fences)} fences, "
            f"{len(self.overpasses)} overpasses.")

        self.active_task = AsfConnectivityTask(
            raster_path=raster.source(),
            origin_xy=origin_xy,
            origin_radius_cells=2,
            fences=fences_geojson,
            overpasses=overpasses_xy,
            window_bounds=(min_x, min_y, max_x, max_y),
            crs_wkt=raster_crs.toWkt(),
            options=options,
        )
        QgsApplication.taskManager().addTask(self.active_task)

    # -----------------------------------------------------------------
    def _show_origin_preview(self, origin_pt, crs):
        """Drop a small marker at the outbreak point so the user can see it."""
        authid = crs.authid() or "EPSG:2056"
        lay = QgsVectorLayer(
            f"Point?crs={authid}", "Outbreak origin", "memory")
        prov = lay.dataProvider()
        f = QgsFeature()
        f.setGeometry(QgsGeometry.fromPointXY(origin_pt))
        prov.addFeature(f)
        lay.updateExtents()
        sym = lay.renderer().symbol()
        sym.setColor(QColor(231, 76, 60))
        sym.setSize(5.0)
        add_wildboar_layer(lay, ZOrder.CENTRE)

    def _collect_habitat_polygons(self, tr_to_raster):
        """Read the user-picked habitat layer and return a list of
        GeoJSON-like polygon dicts in the raster CRS. Returns None if
        no habitat layer is selected (the task then falls back to its
        auto-detected clusters).
        """
        layer = self.dlg.cmbHabitats.currentLayer()
        if layer is None or not isinstance(layer, QgsVectorLayer):
            return None
        if layer.geometryType() != QgsWkbTypes.PolygonGeometry:
            return None

        out = []
        for feat in layer.getFeatures():
            g = feat.geometry()
            if g is None or g.isEmpty():
                continue
            # Transform a copy into the raster CRS before serialising.
            gg = QgsGeometry(g)
            gg.transform(tr_to_raster)
            try:
                out.append(json.loads(gg.asJson()))
            except Exception:
                continue
        QgsMessageLog.logMessage(
            f"Habitats selected: {len(out)} polygons "
            f"from '{layer.name()}'.",
            LOG_TAG, Qgis.Info)
        return out if out else None

    def _warn(self, msg):
        self.dlg.lblStatus.setText(msg)
        QgsMessageLog.logMessage(msg, LOG_TAG, Qgis.Warning)
        QMessageBox.warning(self.dlg, "Wildboar Connectivity", msg)
