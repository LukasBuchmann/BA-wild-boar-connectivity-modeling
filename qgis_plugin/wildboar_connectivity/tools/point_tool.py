# -*- coding: utf-8 -*-
"""
PointSelectionTool - helper map tool for picking a single QgsPointXY on the
QGIS map canvas.

The tool wraps QgsMapToolEmitPoint and adds:

    * A `point_selected(QgsPointXY, str)` signal carrying both the picked
      coordinate and a free-form `role` string (e.g. "start" or "end") so a
      single tool class can serve multiple buttons.
    * An on-canvas QgsVertexMarker that visually anchors the chosen point
      until the next pick or until clear_marker() is called.
    * Automatic restoration of the previous map tool once a point is captured,
      so the user is not stuck in "click to capture" mode.
"""

from qgis.PyQt.QtCore import pyqtSignal, Qt
from qgis.PyQt.QtGui import QColor

from qgis.gui import QgsMapToolEmitPoint, QgsVertexMarker
from qgis.core import QgsPointXY


class PointSelectionTool(QgsMapToolEmitPoint):
    """Map tool that emits a single point click and tags it with a role."""

    # role is "start" or "end" - lets one slot route both signals.
    point_selected = pyqtSignal(QgsPointXY, str)

    def __init__(self, canvas, role: str = "start", color: str = "#2ecc71"):
        super().__init__(canvas)
        self.canvas = canvas
        self.role = role
        self._previous_tool = None
        self._marker = QgsVertexMarker(canvas)
        self._marker.setColor(QColor(color))
        self._marker.setIconSize(12)
        self._marker.setIconType(QgsVertexMarker.ICON_CROSS)
        self._marker.setPenWidth(3)
        self._marker.hide()

    # ------------------------------------------------------------------
    # Activation helpers
    # ------------------------------------------------------------------
    def activate_for_role(self, role: str):
        """Switch the tool into a capture session for the given role."""
        self.role = role
        self._previous_tool = self.canvas.mapTool()
        self.canvas.setMapTool(self)

    def set_point(self, point: QgsPointXY):
        """Programmatically place the marker (used when restoring state)."""
        self._marker.setCenter(point)
        self._marker.show()

    def clear_marker(self):
        self._marker.hide()

    # ------------------------------------------------------------------
    # QgsMapToolEmitPoint overrides
    # ------------------------------------------------------------------
    def canvasReleaseEvent(self, event):
        point = self.toMapCoordinates(event.pos())
        self._marker.setCenter(point)
        self._marker.show()
        self.point_selected.emit(QgsPointXY(point), self.role)

        # Hand the canvas back to whatever tool was active before us
        # so the user is not stuck in capture mode after one click.
        if self._previous_tool is not None:
            self.canvas.setMapTool(self._previous_tool)
            self._previous_tool = None

    def deactivate(self):
        # Explicit override so QGIS always sees this tool as overriding
        # deactivate(), even though the base implementation is sufficient
        # here. Keeps the door open for cleanup logic later without having
        # to remember to add the override first.
        super().deactivate()

    # ------------------------------------------------------------------
    # Clean up the vertex marker when the tool is destroyed
    # ------------------------------------------------------------------
    def __del__(self):
        try:
            if self._marker is not None:
                self.canvas.scene().removeItem(self._marker)
        except Exception:
            # During QGIS teardown the scene may already be gone; ignore.
            pass
