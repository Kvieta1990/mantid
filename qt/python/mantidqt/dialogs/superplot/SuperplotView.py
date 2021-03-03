# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

from qtpy.QtWidgets import QDialog
from qtpy.QtCore import *
from qtpy import uic

import os


class SuperplotView(QDialog):

    _presenter = None

    def __init__(self, presenter):
        self._presenter = presenter
        super().__init__(None)
        self.here = os.path.dirname(os.path.realpath(__file__))
        uic.loadUi(os.path.join(self.here, 'ui/superplot.ui'), self)
