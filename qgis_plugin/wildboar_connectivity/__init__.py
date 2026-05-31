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


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load WildboarConnectivity from main_plugin.py."""
    from .main_plugin import WildboarConnectivity
    return WildboarConnectivity(iface)
