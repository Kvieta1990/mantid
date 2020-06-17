// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#pragma once

#include "DllOption.h"
#include "MantidQtWidgets/Common/IImageInfoWidget.h"
#include "MantidQtWidgets/Common/ImageInfoModel.h"

#include <QTableWidget>

namespace MantidQt {
namespace MantidWidgets {

/**
 * A table widget containing information about the pixel the mouse is over in
 * an image
 */
class EXPORT_OPT_MANTIDQT_COMMON ImageInfoPresenter {

public:
  ImageInfoPresenter(IImageInfoWidget *view);

  std::vector<QString> getInfoList(const double x = DBL_MAX,
                                   const double y = DBL_MAX,
                                   const double z = DBL_MAX);

  void createImageInfoModel(const std::string &ws_type);
  void setWorkspace(const Mantid::API::Workspace_sptr &ws);

  std::shared_ptr<ImageInfoModel> getModel() { return m_model; }

private:
  std::shared_ptr<ImageInfoModel> m_model;
  IImageInfoWidget *m_view;
};

} // namespace MantidWidgets
} // namespace MantidQt
