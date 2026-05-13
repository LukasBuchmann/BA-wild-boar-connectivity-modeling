# -*- coding: utf-8 -*-
"""
WildboarConnectivity v2 - QGIS plugin entry point.

Calculates connectivity between a user-defined start and end point across a
resistance surface using three complementary methods:

    1. Least-Cost Path (LCP)               - single optimal corridor
    2. Circuit Theory (Current Flow)       - pinchpoint / bottleneck map
    3. Graph-based Connectivity            - network visualisation

Begin     : 2026-05-13
Author    : Lukas Buchmann
Email     : buchmluk@students.zhaw.ch
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load the WildboarConnectivity plugin class from main_plugin.py."""
    from .main_plugin import WildboarConnectivity
    return WildboarConnectivity(iface)
