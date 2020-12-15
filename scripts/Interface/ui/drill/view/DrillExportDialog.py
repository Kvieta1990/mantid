# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

import os

from qtpy.QtWidgets import QDialog, QCheckBox
from qtpy import uic


class DrillExportDialog(QDialog):

    UI_FILENAME = "ui/export.ui"

    """
    Export presenter.
    """
    _presenter = None

    def __init__(self, parent=None):
        """
        Create the export dialog.

        Args:
            parent (QWidget): parent widget
        """
        super().__init__(parent)
        self.here = os.path.dirname(os.path.realpath(__file__))
        uic.loadUi(os.path.join(self.here, self.UI_FILENAME), self)
        self.okButton.clicked.connect(self.accept)
        self.cancelButton.clicked.connect(self.reject)
        self.applyButton.clicked.connect(
                lambda : self.accepted.emit()
                )

    def setPresenter(self, presenter):
        """
        Set the prsesenter. This is needed to keep a valid reference to the
        presenter and avoid destruction by the  garage collector.

        Args:
            presenter (DrillExportPresenter): export presenter
        """
        self._presenter = presenter

    def setAlgorithms(self, algorithms):
        """
        Set the algorithms displayed on the dialog.

        Args:
            algorithms (list(str)): list of algorithms
        """
        for i in range(len(algorithms)):
            widget = QCheckBox(algorithms[i], self)
            self.algoList.insertWidget(i, widget)
