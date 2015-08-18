#ifndef MANTIDQTCUSTOMINTERFACES_ENGGDIFFRACTION_ENGGDIFFRACTIONPRESENTER_H_
#define MANTIDQTCUSTOMINTERFACES_ENGGDIFFRACTION_ENGGDIFFRACTIONPRESENTER_H_

#include "MantidAPI/ITableWorkspace_fwd.h"
#include "MantidAPI/MatrixWorkspace_fwd.h"
#include "MantidKernel/System.h"
#include "MantidQtCustomInterfaces/EnggDiffraction/IEnggDiffractionPresenter.h"
#include "MantidQtCustomInterfaces/EnggDiffraction/IEnggDiffractionView.h"
// #include "MantidQtCustomInterfaces/EnggDiffraction/EnggDiffractionModel.h"

#include <boost/scoped_ptr.hpp>

#include <QObject>

class QThread;

namespace MantidQt {
namespace CustomInterfaces {

/**
Presenter for the Enggineering Diffraction GUI (presenter as in the
MVP Model-View-Presenter pattern). In principle, in a strict MVP
setup, signals from the model should always be handled through this
presenter and never go directly to the view, and viceversa.

Copyright &copy; 2015 ISIS Rutherford Appleton Laboratory, NScD
Oak Ridge National Laboratory & European Spallation Source

This file is part of Mantid.

Mantid is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Mantid is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

File change history is stored at: <https://github.com/mantidproject/mantid>
Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
// needs to be dll-exported for the tests
class DLLExport EnggDiffractionPresenter : public QObject,
                                           public IEnggDiffractionPresenter {
  // Q_OBJECT for 'connect' with thread/worker
  Q_OBJECT

public:
  /// Default constructor - normally used from the concrete view
  EnggDiffractionPresenter(IEnggDiffractionView *view);
  virtual ~EnggDiffractionPresenter();

  virtual void notify(IEnggDiffractionPresenter::Notification notif);

  // this is the hard work that a worker / thread will run
  void doNewCalibration(const std::string &outFilename,
                        const std::string &vanNo, const std::string &ceriaNo);

protected:
  void initialize();

  /// clean shut down of model, view, etc.
  void cleanup();

  void processStart();
  void processLoadExistingCalib();
  void processCalcCalib();
  void processLogMsg();
  void processInstChange();
  void processShutDown();

protected slots:
  void calibrationFinished();

private:
  std::string buildCalibrateSuggestedFilename(const std::string &vanNo,
                                              const std::string &ceriaNo) const;

  virtual void startAsynCalibWorker(const std::string &outFilename,
                                    const std::string &vanNo,
                                    const std::string &ceriaNo);

  void doCalib(const EnggDiffCalibSettings &cs, const std::string &vanNo,
               const std::string &ceriaNo, const std::string &outFilename);

  void parseCalibrateFilename(const std::string &path, std::string &instName,
                              std::string &vanNo, std::string &ceriaNo);

  void
  loadVanadiumPrecalcWorkspaces(const std::string &preIntegFilename,
                                const std::string &preCurvesFilename,
                                Mantid::API::ITableWorkspace_sptr &vanIntegWS,
                                Mantid::API::MatrixWorkspace_sptr &vanCurvesWS);

  void calcVanadiumWorkspaces(const std::string &vanNo,
                              Mantid::API::ITableWorkspace_sptr &vanIntegWS,
                              Mantid::API::MatrixWorkspace_sptr &vanCurvesWS);

  void findPrecalcVanadiumCorrFilenames(const std::string &vanNo,
                                        const std::string &inputDirCalib,
                                        std::string &preIntegFilename,
                                        std::string &preCurvesFilename,
                                        bool &found);

  void loadOrCalcVanadiumWorkspaces(
      const std::string &vanNo, const std::string &inputDirCalib,
      Mantid::API::ITableWorkspace_sptr &vanIntegWS,
      Mantid::API::MatrixWorkspace_sptr &vanCurvesWS, bool forceRecalc);

  /// string to use for ENGINX file names (as a prefix, etc.)
  const static std::string g_enginxStr;

  /// whether to allow users to give the output calibration filename
  const static bool g_askUserCalibFilename;

  QThread *m_workerThread;

  /// true if the last calibration completed successfully
  bool m_calibFinishedOK;

  /// Associated view for this presenter (MVP pattern)
  IEnggDiffractionView *const m_view;

  /// Associated model for this presenter (MVP pattern)
  // const boost::scoped_ptr<EnggDiffractionModel> m_model;
};

} // namespace CustomInterfaces
} // namespace MantidQt

#endif // MANTIDQTCUSTOMINTERFACES_ENGGDIFFRACTION_ENGGDIFFRACTIONPRESENTER_H_
