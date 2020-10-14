// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidQtWidgets/Common/MuonFunctionBrowser.h"
#include "MantidQtWidgets/Common/SelectFunctionDialog.h"

namespace MantidQt {
namespace MantidWidgets {

/**
 * Constructor
 * @param parent :: The parent widget.
 * @param multi  :: Option to use the browser for multi-dataset fitting.
 */
MuonFunctionBrowser::MuonFunctionBrowser(QWidget *parent, bool multi)
    : FunctionBrowser(parent, multi, {"Muon", "General", "Background"}) {}

/**
 * Destructor
 */
MuonFunctionBrowser::~MuonFunctionBrowser() {}

} // namespace MantidWidgets
} // namespace MantidQt
