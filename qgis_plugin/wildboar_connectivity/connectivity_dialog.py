# -*- coding: utf-8 -*-
"""Dialog loader for the Wildboar Connectivity (ASF) plugin.

The dialog's layout lives in connectivity_dialog_base.ui (a Qt Designer
file), not in Python. uic.loadUiType() reads that file at import time and
builds a Python class from it, so the form can be edited visually in Qt
Designer without ever touching this module.
"""

import os

from qgis.PyQt import uic, QtWidgets

# Build the form class from the .ui file that sits next to this module.
FORM_CLASS, _ = uic.loadUiType(
    os.path.join(os.path.dirname(__file__), 'connectivity_dialog_base.ui'))


class WildboarConnectivityDialog(QtWidgets.QDialog, FORM_CLASS):
    """Main plugin dialog. All widgets are defined in the .ui file and
    become attributes of `self` once setupUi() runs."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
