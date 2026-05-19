# -*- coding: utf-8 -*-
"""
FenceDrawingTool - polyline-drawing map tool for placing virtual fences
on the canvas.

Usage:
    * Left-click to add a vertex.
    * Move the mouse to see a rubber-band preview of the next segment.
    * Right-click or double-click to finish the polyline.
    * Press Escape to cancel without emitting.

When a polyline is finished, the tool emits `fence_completed(QgsGeometry)`
carrying the polyline in the canvas CRS. The plugin transforms it into
the raster CRS before storing it.
"""

from qgis.PyQt.QtCore import Qt, pyqtSignal
from qgis.PyQt.QtGui import QColor

from qgis.core import QgsGeometry, QgsPointXY, QgsWkbTypes
from qgis.gui import QgsMapTool, QgsRubberBand


class FenceDrawingTool(QgsMapTool):
    """Polyline-drawing map tool that emits the finished geometry."""

    fence_completed = pyqtSignal(QgsGeometry)

    def __init__(self, canvas, color: str = "#e74c3c"):
        super().__init__(canvas)
        self.canvas = canvas

        self.points = []                 # confirmed vertices
        self._previous_tool = None

        # Rubber band shown while drawing.
        self.rb = QgsRubberBand(canvas, QgsWkbTypes.LineGeometry)
        self.rb.setColor(QColor(color))
        self.rb.setWidth(3)

    # ------------------------------------------------------------------
    def activate_for_drawing(self):
        """Switch the canvas to this tool and remember the previous one."""
        self._previous_tool = self.canvas.mapTool()
        self.canvas.setMapTool(self)

    def reset(self):
        self.points = []
        self.rb.reset(QgsWkbTypes.LineGeometry)

    # ------------------------------------------------------------------
    # QgsMapTool overrides
    # ------------------------------------------------------------------
    def canvasPressEvent(self, event):
        if event.button() == Qt.LeftButton:
            pt = self.toMapCoordinates(event.pos())
            self.points.append(QgsPointXY(pt))
            # First click resets the rubber band; subsequent clicks extend.
            if len(self.points) == 1:
                self.rb.reset(QgsWkbTypes.LineGeometry)
                self.rb.addPoint(pt, True)
            else:
                self.rb.addPoint(pt, True)
        elif event.button() == Qt.RightButton:
            self._finish()

    def canvasMoveEvent(self, event):
        if not self.points:
            return
        # Live preview: the last rubber-band point follows the mouse.
        pt = self.toMapCoordinates(event.pos())
        n = self.rb.numberOfVertices()
        if n > len(self.points):
            self.rb.removeLastPoint()
        self.rb.addPoint(pt, True)

    def canvasDoubleClickEvent(self, event):
        self._finish()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            self.reset()
            self._restore_previous_tool()

    # ------------------------------------------------------------------
    def _finish(self):
        # Trim the floating preview vertex if present.
        if self.rb.numberOfVertices() > len(self.points):
            self.rb.removeLastPoint()
        if len(self.points) >= 2:
            geom = QgsGeometry.fromPolylineXY(list(self.points))
            if not geom.isEmpty():
                self.fence_completed.emit(geom)
        self.reset()
        self._restore_previous_tool()

    def _restore_previous_tool(self):
        if self._previous_tool is not None:
            self.canvas.setMapTool(self._previous_tool)
            self._previous_tool = None

    # ------------------------------------------------------------------
    def __del__(self):
        try:
            self.canvas.scene().removeItem(self.rb)
        except Exception:
            pass
