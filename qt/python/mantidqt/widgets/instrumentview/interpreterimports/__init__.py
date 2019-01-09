# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
from __future__ import (absolute_import, division, print_function, unicode_literals)

# local imports
from mantidqt.utils.qt import import_qt

# Import widget from C++ wrappers
InstrumentWidgetEncoder = import_qt('._instrumentview', 'mantidqt.widgets.instrumentview',
                                    'InstrumentWidgetEncoder')
InstrumentWidgetDecoder = import_qt('._instrumentview', 'mantidqt.widgets.instrumentview',
                                    'InstrumentWidgetDecoder')
