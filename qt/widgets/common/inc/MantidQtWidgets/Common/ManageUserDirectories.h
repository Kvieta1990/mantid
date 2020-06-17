// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "DllOption.h"
#include "MantidQtWidgets/Common/MantidDialog.h"
#include "ui_ManageUserDirectories.h"
#include <QDialog>

namespace MantidQt {
namespace API {

class EXPORT_OPT_MANTIDQT_COMMON ManageUserDirectories
    : public MantidQt::API::MantidDialog {
  Q_OBJECT

public:
  ManageUserDirectories(QWidget *parent = nullptr);
  ~ManageUserDirectories() override;
  static void openManageUserDirectories();
  void closeEvent(QCloseEvent *event) override;

private:
  virtual void initLayout();
  void loadProperties();
  void saveProperties();
  void appendSlashIfNone(QString &path) const;
  QListWidget *listWidget(QObject *object);

private slots:
  void helpClicked();
  void cancelClicked();
  void confirmClicked();
  void addDirectory();
  void browseToDirectory();
  void remDir();
  void moveUp();
  void moveDown();
  void selectSaveDir();

private:
  Ui::ManageUserDirectories m_uiForm;
  QString m_userPropFile;
};

} // namespace API
} // namespace MantidQt
