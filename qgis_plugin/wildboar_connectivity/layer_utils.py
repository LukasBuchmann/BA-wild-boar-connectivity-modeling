# -*- coding: utf-8 -*-
"""
Shared layer-management helpers for the Wildboar Connectivity plugin.

Each output layer is added to the project with `addToLegend=False` (so
QGIS doesn't auto-create a tree node) and then inserted into the layer
tree at the correct z-position via `QgsLayerTreeGroup.insertLayer()`.

Canonical top-to-bottom layer order (lower z = higher in the tree =
rendered on top):

    0  Outbreak origin point
    1  Cluster destination dots
    2  Fences / overpasses preview
    3  Dispersal corridors (LCP lines)
    4  Infection risk raster
    5  Dispersal corridor probability raster (Circuit theory)
    6  Random-walk density raster
    7  Contamination probability raster (iSSF-IBMM)
    8  Swissimage basemap (60 % opacity)

Untagged layers (resistance raster, habitat polygons, etc.) stay where
they are, below the tagged stack.
"""

from qgis.core import (
    QgsProject,
    QgsRasterLayer,
)


WB_ZORDER_KEY  = "wildboar_z_order"
WB_BASEMAP_TAG = "wildboar_basemap"


class ZOrder:
    """Lower z = higher in the QGIS layer tree = rendered on top."""
    CENTRE        = 0   # Outbreak origin marker
    ANCHORS       = 1   # Cluster destination dots
    FENCES        = 2   # Fences / overpasses interactive preview
    LCP_TRAFFIC   = 3   # Dispersal corridor lines (LCPs)
    RISK          = 4   # Infection risk raster
    PINCHPOINTS   = 5   # Circuit-theory current (pinchpoints) raster
    WALK_DENSITY  = 6   # BCRW random-walk density raster
    IBMM_DENSITY  = 7   # iSSF-IBMM contamination probability raster
    BASEMAP       = 8   # Swissimage basemap — always at the bottom


def _correct_index_for_z(z: int) -> int:
    """Return the tree-child index where a new layer with z should be inserted.

    Walks root.children() top-to-bottom and returns the index of the first
    child that should sit below the new layer (higher z or untagged).
    """
    root = QgsProject.instance().layerTreeRoot()
    for i, child in enumerate(root.children()):
        layer = getattr(child, "layer", lambda: None)()
        if layer is None:
            return i
        other_z = layer.customProperty(WB_ZORDER_KEY)
        if other_z in (None, ""):
            return i
        try:
            if z < int(other_z):
                return i
        except (TypeError, ValueError):
            return i
    return len(root.children())


def add_wildboar_layer(layer, z):
    """Add *layer* to the project and insert it at the correct z-position."""
    from qgis.core import Qgis, QgsMessageLog
    project = QgsProject.instance()
    layer.setCustomProperty(WB_ZORDER_KEY, int(z))
    project.addMapLayer(layer, addToLegend=False)
    project.layerTreeRoot().insertLayer(_correct_index_for_z(int(z)), layer)
    QgsMessageLog.logMessage(
        f"add_wildboar_layer: '{layer.name()}' z={z} valid={layer.isValid()}",
        "wildboar-asf", Qgis.Info)
    return layer


def ensure_swissimage(opacity: float = 0.6):
    """Ensure the swisstopo Swissimage WMTS basemap is present in the project."""
    project = QgsProject.instance()
    root = project.layerTreeRoot()

    for layer in project.mapLayers().values():
        if layer.customProperty(WB_BASEMAP_TAG) == "1":
            try:
                layer.setOpacity(opacity)
            except Exception:
                pass
            layer.setCustomProperty(WB_ZORDER_KEY, int(ZOrder.BASEMAP))
            if root.findLayer(layer.id()) is None:
                root.insertLayer(_correct_index_for_z(int(ZOrder.BASEMAP)), layer)
            return layer

    uri = (
        "type=xyz"
        "&url=https://wmts.geo.admin.ch/1.0.0/"
        "ch.swisstopo.swissimage/default/current/3857/{z}/{x}/{y}.jpeg"
        "&zmax=18&zmin=0"
    )
    layer = QgsRasterLayer(uri, "Swissimage (swisstopo)", "wms")
    if not layer.isValid():
        return None
    layer.setCustomProperty(WB_BASEMAP_TAG, "1")
    try:
        layer.setOpacity(opacity)
    except Exception:
        pass
    add_wildboar_layer(layer, ZOrder.BASEMAP)
    return layer


def remove_all_wildboar_layers(keep_basemap: bool = False):
    """Remove every plugin-tagged layer from the project.

    Parameters
    ----------
    keep_basemap : bool
        If True, the Swissimage basemap layer is preserved.
    """
    from qgis.core import Qgis, QgsMessageLog
    project = QgsProject.instance()

    ids_to_remove = []
    for lid, lyr in project.mapLayers().items():
        if lyr.customProperty(WB_ZORDER_KEY) in (None, ""):
            continue
        if keep_basemap and lyr.customProperty(WB_BASEMAP_TAG) == "1":
            continue
        ids_to_remove.append(lid)

    if ids_to_remove:
        project.removeMapLayers(ids_to_remove)
        QgsMessageLog.logMessage(
            f"Removed {len(ids_to_remove)} wildboar layer(s) "
            f"(keep_basemap={keep_basemap})",
            "wildboar-asf", Qgis.Info)
