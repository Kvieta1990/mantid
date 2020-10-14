// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2014 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "ALFCustomInstrumentModel.h"
#include "ALFCustomInstrumentPresenter.h"
#include "ALFCustomInstrumentView.h"
#include "DllConfig.h"
#include "MantidQtWidgets/Common/ObserverPattern.h"
#include "MantidQtWidgets/Common/UserSubWindow.h"
#include "MantidQtWidgets/InstrumentView/PlotFitAnalysisPanePresenter.h"

namespace MantidQt {
namespace CustomInterfaces {
/** ALFView : Custom interface for looking at ALF data
 */
class MANTIDQT_DIRECT_DLL ALFView : public API::UserSubWindow {
  Q_OBJECT

public:
  ALFView(QWidget *parent = nullptr);
  ~ALFView() { delete m_presenter; };
  static std::string name() { return "ALF View"; }
  static QString categoryInfo() { return "Direct"; }

protected:
  void initLayout() override;

private:
  ALFCustomInstrumentView *m_view;
  ALFCustomInstrumentModel *m_model;
  ALFCustomInstrumentPresenter *m_presenter;
  MantidWidgets::PlotFitAnalysisPanePresenter *m_analysisPane;
};
} // namespace CustomInterfaces
} // namespace MantidQt
