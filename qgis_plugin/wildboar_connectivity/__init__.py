# -*- coding: utf-8 -*-
"""
WildboarConnectivity (ASF) - QGIS plugin entry point.

Habitat-free decision-support tool for African Swine Fever spread:

    1. Load a resistance raster (Keeley-transformed iSSF surface).
    2. Click the outbreak origin on the map.
    3. Optionally draw fences and/or place wildlife overpasses.
    4. Press Run. Outputs:
         - Infection-risk raster (resistant-kernel, green/yellow/red)
         - LCP corridors (least-cost paths to nearest low-resistance clusters)
         - Pinchpoints raster (Circuit theory, Circuitscape.jl or scipy)
         - BCRW random-walk density (optional)
         - iSSF-IBMM contamination probability (optional)

Author : Lukas Buchmann <buchmluk@students.zhaw.ch>
"""


# QGIS calls this function by name when it loads the plugin. The name and
# the camelCase signature are fixed by the QGIS plugin API, which is why the
# two comments below silence the linter warnings about the naming style.
# The import is done inside the function (instead of at the top of the file)
# so that QGIS can discover the plugin without first importing the heavy
# main_plugin module, which in turn pulls in all the analysis dependencies.
# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Create and return the plugin instance, handing it the QGIS interface."""
    from .main_plugin import WildboarConnectivity
    return WildboarConnectivity(iface)
