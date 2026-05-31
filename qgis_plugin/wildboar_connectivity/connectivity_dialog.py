# -*- coding: utf-8 -*-
"""Dialog loader for the Wildboar Connectivity (ASF) plugin."""

import os

from qgis.PyQt import uic, QtWidgets

FORM_CLASS, _ = uic.loadUiType(
    os.path.join(os.path.dirname(__file__), 'connectivity_dialog_base.ui'))


class WildboarConnectivityDialog(QtWidgets.QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
