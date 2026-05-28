# -*- coding: utf-8 -*-
"""
Shared layer-management helpers for the Wildboar Connectivity plugin.

Each output layer is added to the project with `addToLegend=False` (so
QGIS doesn't auto-create a tree node) and then inserted into the layer
tree at the correct z-position via `QgsLayerTreeGroup.insertLayer()`.

That single-step add avoids the disappearing-layer bug that earlier
versions hit: the layer-tree <-> registry bridge was auto-removing
layers from the project whenever a tree node was removed (and our
old reorder pass removed every tagged node before re-inserting them).
Inserting at the right index from the start means we never detach a
node, so the bridge has nothing to react to.

Canonical top-to-bottom layer order (lower z = higher in the tree):

    0  Google Satellite Hybrid (30 % opacity)
    1  Selected habitats
    2  LCP traffic raster
    3  Pinchpoints raster
    4  Habitat graph - nodes
    5  Habitat graph - edges
    6  Analysis centre (point)
    7  Analysis range (buffer)

Untagged layers (resistance raster, full core-habitat input, etc.) stay
where they are, below the tagged stack.
"""

from qgis.core import (
    QgsProject,
    QgsRasterLayer,
)


# Custom-property keys
WB_ZORDER_KEY = "wildboar_z_order"
WB_GOOGLE_TAG = "wildboar_google_satellite"


class ZOrder:
    """Lower number = higher in the tree = drawn on top."""
    GOOGLE             = 0
    SELECTED_HABITATS  = 1
    LCP_TRAFFIC        = 2
    PINCHPOINTS        = 3
    IBMM_DENSITY       = 4   # iSSF-IBMM agent visit density
    GRAPH_NODES        = 5
    GRAPH_EDGES        = 6
    CENTRE             = 7
    RANGE              = 8


def _correct_index_for_z(z: int) -> int:
    """Compute the tree-child index where a new layer with this z should
    go. Walks root.children() top-to-bottom and returns the index of the
    first child that should be below us:

        - an existing tagged layer with a higher z value, OR
        - the first untagged layer (tagged layers always sit above
          untagged layers in our stack).
    """
    root = QgsProject.instance().layerTreeRoot()
    children = root.children()
    for i, child in enumerate(children):
        layer = getattr(child, "layer", lambda: None)()
        if layer is None:
            # A group node, etc. - treat as untagged.
            return i
        other_z = layer.customProperty(WB_ZORDER_KEY)
        if other_z in (None, ""):
            return i                       # first untagged: we go above it
        try:
            other_z_int = int(other_z)
        except (TypeError, ValueError):
            return i
        if z < other_z_int:
            return i                       # lower z = above this child
        # else: continue, we go further down
    return len(children)                   # append at the bottom


def add_wildboar_layer(layer, z):
    """Add a layer to the project AND to the tree at the right position.

    The layer becomes visible immediately - region-preview layers
    (centre, range, selected habitats, Google Satellite) appear the
    moment the user clicks Run, even while the analysis task is still
    crunching in the background.
    """
    from qgis.core import Qgis, QgsMessageLog
    project = QgsProject.instance()
    layer.setCustomProperty(WB_ZORDER_KEY, int(z))
    # addToLegend=False: don't let QGIS pick a tree position for us.
    project.addMapLayer(layer, addToLegend=False)

    target_idx = _correct_index_for_z(int(z))
    project.layerTreeRoot().insertLayer(target_idx, layer)

    QgsMessageLog.logMessage(
        f"add_wildboar_layer: '{layer.name()}' z={z} idx={target_idx} "
        f"valid={layer.isValid()}",
        "wildboar-v2", Qgis.Info)
    return layer


def ensure_google_satellite(opacity: float = 0.3):
    """Make sure a Swissimage WMTS basemap is present in the project.

    Uses the official swisstopo Swissimage tile service - the standard
    Swiss orthoimagery, far more appropriate for a Swiss thesis than
    a Google product. The function name and the tag string
    `WB_GOOGLE_TAG` are kept for backward compatibility so layers
    added by older plugin versions are still recognised and reused
    rather than duplicated.
    """
    project = QgsProject.instance()
    root = project.layerTreeRoot()

    for layer in project.mapLayers().values():
        if layer.customProperty(WB_GOOGLE_TAG) == "1":
            try:
                layer.setOpacity(opacity)
            except Exception:
                pass
            layer.setCustomProperty(WB_ZORDER_KEY, int(ZOrder.GOOGLE))
            if root.findLayer(layer.id()) is None:
                target_idx = _correct_index_for_z(int(ZOrder.GOOGLE))
                root.insertLayer(target_idx, layer)
            return layer

    # Swissimage WMTS endpoint (Web Mercator tiles, public, free). The
    # URL has no query parameters, so no percent-encoding is needed -
    # the `{z}/{x}/{y}` placeholders are stamped in by QGIS at draw time.
    url_template = (
        "https://wmts.geo.admin.ch/1.0.0/"
        "ch.swisstopo.swissimage/default/current/3857/"
        "{z}/{x}/{y}.jpeg"
    )
    uri = f"type=xyz&url={url_template}&zmax=18&zmin=0"

    layer = QgsRasterLayer(uri, "Swissimage (swisstopo)", "wms")
    if not layer.isValid():
        return None
    layer.setCustomProperty(WB_GOOGLE_TAG, "1")
    try:
        layer.setOpacity(opacity)
    except Exception:
        pass
    add_wildboar_layer(layer, ZOrder.GOOGLE)
    return layer


def reorder_wildboar_layers():
    """Back-compat no-op.

    Each layer is now inserted at the correct z-position when it is
    added (see add_wildboar_layer), so a separate reorder pass is no
    longer needed - and removing one would re-trigger the bridge bug
    that made the layers disappear.
    """
    from qgis.core import Qgis, QgsMessageLog
    project = QgsProject.instance()
    final_count = sum(
        1 for layer in project.mapLayers().values()
        if layer.customProperty(WB_ZORDER_KEY) not in (None, "")
    )
    QgsMessageLog.logMessage(
        f"reorder (no-op): {final_count} wildboar layers in project",
        "wildboar-v2", Qgis.Info)


def remove_all_wildboar_layers(keep_google: bool = False):
    """Delete every layer created by this plugin from the current project.

    A layer is "ours" if it carries the `wildboar_z_order` custom
    property (set by `add_wildboar_layer`). User-loaded layers
    (resistance raster, core-habitat polygons, anything else they
    opened) are never touched.

    Parameters
    ----------
    keep_google : bool
        If True, the Google Satellite Hybrid tile layer is preserved.
        Useful for the pre-run cleanup (each run swaps in fresh outputs
        but the user almost certainly wants to keep the basemap).
    """
    from qgis.core import Qgis, QgsMessageLog
    project = QgsProject.instance()

    ids_to_remove = []
    for lid, lyr in project.mapLayers().items():
        z = lyr.customProperty(WB_ZORDER_KEY)
        if z in (None, ""):
            continue
        if keep_google and lyr.customProperty(WB_GOOGLE_TAG) == "1":
            continue
        ids_to_remove.append(lid)

    if ids_to_remove:
        project.removeMapLayers(ids_to_remove)
        QgsMessageLog.logMessage(
            f"removed {len(ids_to_remove)} wildboar layer(s) "
            f"(keep_google={keep_google})",
            "wildboar-v2", Qgis.Info)
