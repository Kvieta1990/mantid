#include "MantidAPI/ITableWorkspace.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidQtAPI/AlgorithmInputHistory.h"
#include "MantidQtAPI/AlgorithmRunner.h"
#include "MantidQtAPI/PythonRunner.h"
// #include "MantidQtCustomInterfaces/EnggDiffraction/EnggDiffractionModel.h"
#include "MantidQtCustomInterfaces/EnggDiffraction/EnggDiffractionPresenter.h"
#include "MantidQtCustomInterfaces/EnggDiffraction/EnggDiffractionPresWorker.h"
#include "MantidQtCustomInterfaces/EnggDiffraction/IEnggDiffractionView.h"

#include <boost/lexical_cast.hpp>

#include <Poco/File.h>
#include <Poco/Path.h>

#include <QThread>

using namespace Mantid::API;
using namespace MantidQt::CustomInterfaces;

namespace MantidQt {
namespace CustomInterfaces {

namespace {
Mantid::Kernel::Logger g_log("EngineeringDiffractionGUI");
}

const std::string EnggDiffractionPresenter::g_enginxStr = "ENGINX";
// discouraged at the moment
const bool EnggDiffractionPresenter::g_askUserCalibFilename = false;

EnggDiffractionPresenter::EnggDiffractionPresenter(IEnggDiffractionView *view)
    : m_workerThread(NULL), m_calibFinishedOK(false),
      m_view(view) /*, m_model(new EnggDiffractionModel()), */ {
  if (!m_view) {
    throw std::runtime_error(
        "Severe inconsistency found. Presenter created "
        "with an empty/null view (engineeering diffraction interface). "
        "Cannot continue.");
  }
}

EnggDiffractionPresenter::~EnggDiffractionPresenter() { cleanup(); }

/**
 * Close open sessions, kill threads etc., save settings, etc. for a
 * graceful window close/destruction
 */
void EnggDiffractionPresenter::cleanup() {
  // m_model->cleanup();

  // this may still be running
  if (m_workerThread) {
    if (m_workerThread->isRunning()) {
      g_log.notice() << "A calibration process is currently running, shutting "
                        "it down immediately..." << std::endl;
      m_workerThread->wait(10);
    }
    delete m_workerThread;
    m_workerThread = NULL;
  }
}

void EnggDiffractionPresenter::notify(
    IEnggDiffractionPresenter::Notification notif) {

  switch (notif) {

  case IEnggDiffractionPresenter::Start:
    processStart();
    break;

  case IEnggDiffractionPresenter::LoadExistingCalib:
    processLoadExistingCalib();
    break;

  case IEnggDiffractionPresenter::CalcCalib:
    processCalcCalib();
    break;

  case IEnggDiffractionPresenter::LogMsg:
    processLogMsg();
    break;

  case IEnggDiffractionPresenter::InstrumentChange:
    processInstChange();

  case IEnggDiffractionPresenter::ShutDown:
    processShutDown();
    break;
  }
}

void EnggDiffractionPresenter::processStart() {
  EnggDiffCalibSettings cs = m_view->currentCalibSettings();
}

void EnggDiffractionPresenter::processLoadExistingCalib() {
  EnggDiffCalibSettings cs = m_view->currentCalibSettings();

  std::string fname = m_view->askExistingCalibFilename();
  if (fname.empty()) {
    return;
  }

  std::string instName, vanNo, ceriaNo;
  try {
    parseCalibrateFilename(fname, instName, vanNo, ceriaNo);
  } catch (std::invalid_argument &ia) {
    m_view->userWarning("Invalid calibration filename : " + fname, ia.what());
    return;
  }

  m_view->newCalibLoaded(vanNo, ceriaNo, fname);
}

void EnggDiffractionPresenter::processCalcCalib() {
  std::string vanNo = m_view->newVanadiumNo();
  if (vanNo.empty()) {
    m_view->userWarning(
        "Error in the inputs",
        "The Vanadium number cannot be empty and must be an integer number.");
    return;
  }
  std::string ceriaNo = m_view->newCeriaNo();
  if (ceriaNo.empty()) {
    m_view->userWarning(
        "Error in the inputs",
        "The Ceria number cannot be empty and must be an integer number.");
    return;
  }

  EnggDiffCalibSettings cs = m_view->currentCalibSettings();
  const std::string pixelCalib = cs.m_pixelCalibFilename;
  if (pixelCalib.empty()) {
    m_view->userWarning(
        "Error in the inputs",
        "You need to set a pixel (full) calibration in settings.");
    return;
  }
  const std::string templGSAS = cs.m_templateGSAS_PRM;
  if (templGSAS.empty()) {
    m_view->userWarning(
        "Error in the inputs",
        "You need to set a template calibration file for GSAS in settings.");
    return;
  }

  g_log.notice() << "EnggDiffraction GUI: starting new calibration. This may "
                    "take a few seconds... " << std::endl;

  const std::string sugg = buildCalibrateSuggestedFilename(vanNo, ceriaNo);
  std::string outFilename;
  if (!g_askUserCalibFilename) {
    outFilename = sugg;
  } else {
    outFilename = m_view->askNewCalibrationFilename(sugg);
    if (outFilename.empty()) {
      return;
    }
    std::string instName = "";

    // make sure it follows the rules
    try {
      parseCalibrateFilename(outFilename, instName, vanNo, ceriaNo);
    } catch (std::invalid_argument &ia) {
      m_view->userWarning("Invalid output calibration filename: " + outFilename,
                          ia.what());
      return;
    }
  }

  m_view->enableCalibrateActions(false);
  // alternatively, this would be GUI-blocking:
  // doNewCalibration(outFilename, vanNo, ceriaNo);
  // calibrationFinished()
  startAsynCalibWorker(outFilename, vanNo, ceriaNo);
}

/**
 * Start the calibration work without blocking the GUI
 *
 * @param outFilename name for the output GSAS calibration file
 * @param vanNo vanadium run number
 * @param ceriaNo ceria run number
 */
void
EnggDiffractionPresenter::startAsynCalibWorker(const std::string &outFilename,
                                               const std::string &vanNo,
                                               const std::string &ceriaNo) {
  delete m_workerThread;
  m_workerThread = new QThread(this);
  EnggDiffWorker *worker =
      new EnggDiffWorker(this, outFilename, vanNo, ceriaNo);
  worker->moveToThread(m_workerThread);

  connect(m_workerThread, SIGNAL(started()), worker, SLOT(calibrate()));
  connect(worker, SIGNAL(finished()), this, SLOT(calibrationFinished()));
  // early delete of thread and worker
  connect(m_workerThread, SIGNAL(finished()), m_workerThread,
          SLOT(deleteLater()), Qt::DirectConnection);
  connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
  m_workerThread->start();
}

void EnggDiffractionPresenter::calibrationFinished() {
  if (!m_view)
    return;

  m_view->enableCalibrateActions(true);
  if (!m_calibFinishedOK) {
    g_log.warning() << "The cablibration did not finished correctly."
                    << std::endl;
  } else {
    const std::string vanNo = m_view->newVanadiumNo();
    const std::string ceriaNo = m_view->newCeriaNo();
    const std::string outFilename =
        buildCalibrateSuggestedFilename(vanNo, ceriaNo);
    m_view->newCalibLoaded(vanNo, ceriaNo, outFilename);
    g_log.notice()
        << "Cablibration finished and ready as 'current calibration'."
        << std::endl;
  }
  if (m_workerThread) {
    delete m_workerThread;
    m_workerThread = NULL;
  }
}

void EnggDiffractionPresenter::doNewCalibration(const std::string &outFilename,
                                                const std::string &vanNo,
                                                const std::string &ceriaNo) {
  g_log.notice() << "Generating new calibration file: " << outFilename
                 << std::endl;

  EnggDiffCalibSettings cs = m_view->currentCalibSettings();
  Mantid::Kernel::ConfigServiceImpl &conf =
      Mantid::Kernel::ConfigService::Instance();
  const std::vector<std::string> tmpDirs = conf.getDataSearchDirs();
  // in principle, the run files will be found from 'DirRaw', and the
  // pre-calculated Vanadium corrections from 'DirCalib'
  if (!cs.m_inputDirCalib.empty() && Poco::File(cs.m_inputDirCalib).exists()) {
    conf.appendDataSearchDir(cs.m_inputDirCalib);
  }
  if (!cs.m_inputDirRaw.empty() && Poco::File(cs.m_inputDirRaw).exists())
    conf.appendDataSearchDir(cs.m_inputDirRaw);

  try {
    doCalib(cs, vanNo, ceriaNo, outFilename);
    m_calibFinishedOK = true;
  } catch (std::runtime_error &) {
    g_log.error() << "The calibration calculations failed. One of the "
                     "algorithms did not execute correctly. See log messages "
                     "for details. " << std::endl;
  } catch (std::invalid_argument &) {
    g_log.error()
        << "The calibration calculations failed. Some input properties "
           "were not valid. See log messages for details. " << std::endl;
  }
  // restore normal data search paths
  conf.setDataSearchDirs(tmpDirs);
}

void EnggDiffractionPresenter::processLogMsg() {
  std::vector<std::string> msgs = m_view->logMsgs();
  for (size_t i = 0; i < msgs.size(); i++) {
    g_log.information() << msgs[i] << std::endl;
  }
}

void EnggDiffractionPresenter::processInstChange() {
  const std::string err = "Changing instrument is not supported!";
  g_log.error() << err << std::endl;
  m_view->userError("Fatal error", err);
}

void EnggDiffractionPresenter::processShutDown() {
  m_view->saveSettings();
  cleanup();
}

/**
 * Calculate a calibration, responding the the "new calibration"
 * action/button.
 *
 * @param cs user settings
 * @param vanNo Vanadium run number
 * @param ceriaNo Ceria run number
 * @param outFilename output filename chosen by the user
 */
void EnggDiffractionPresenter::doCalib(const EnggDiffCalibSettings &cs,
                                       const std::string &vanNo,
                                       const std::string &ceriaNo,
                                       const std::string &outFilename) {
  ITableWorkspace_sptr vanIntegWS;
  MatrixWorkspace_sptr vanCurvesWS;
  MatrixWorkspace_sptr ceriaWS;

  loadOrCalcVanadiumWorkspaces(vanNo, cs.m_inputDirCalib, vanIntegWS,
                               vanCurvesWS, cs.m_forceRecalcOverwrite);

  const std::string instStr = m_view->currentInstrument();
  try {
    auto load = Algorithm::fromString("Load");
    load->initialize();
    load->setPropertyValue("Filename", instStr + ceriaNo);
    std::string ceriaWSName = "engggui_calibration_sample_ws";
    load->setPropertyValue("OutputWorkspace", ceriaWSName);
    load->execute();

    AnalysisDataServiceImpl &ADS = Mantid::API::AnalysisDataService::Instance();
    ceriaWS = ADS.retrieveWS<MatrixWorkspace>(ceriaWSName);
  } catch (std::runtime_error &re) {
    g_log.error()
        << "Error while loading calibration sample data. "
           "Could not run the algorithm Load succesfully for the calibration "
           "sample (run number: " +
               ceriaNo + "). Error description: " + re.what() +
               " Please check also the previous log messages for details.";
    throw;
  }

  // Bank 1 and 2 - ENGIN-X
  const size_t numBanks = 2;
  std::vector<double> difc, tzero;
  difc.resize(numBanks);
  tzero.resize(numBanks);
  for (size_t i = 0; i < difc.size(); i++) {
    auto alg = Algorithm::fromString("EnggCalibrate");
    try {
      alg->initialize();
      alg->setProperty("InputWorkspace", ceriaWS);
      alg->setProperty("VanIntegrationWorkspace", vanIntegWS);
      alg->setProperty("VanCurvesWorkspace", vanCurvesWS);
      alg->setPropertyValue("Bank", boost::lexical_cast<std::string>(i + 1));
      // TODO: figure out what should be done about the list of expected peaks
      // to EnggCalibrate => it should be a default, as in EnggFitPeaks, that
      // should be fixed in a nother ticket/issue
      alg->setPropertyValue(
          "ExpectedPeaks",
          "3.1243, 2.7057, 1.9132, 1.6316, 1.5621, "
          "1.3529, 1.2415, 1.2100, 1.1046, 1.0414, 0.9566, 0.9147, 0.9019, "
          "0.8556, 0.8252, 0.8158, 0.7811");
      alg->setPropertyValue("OutputParametersTableName",
                            "engggui_calibration_bank_" +
                                boost::lexical_cast<std::string>(i + 1));
      alg->execute();
    } catch (std::runtime_error &re) {
      g_log.error() << "Error in calibration. ",
          "Could not run the algorithm EnggCalibrate succesfully for bank " +
              boost::lexical_cast<std::string>(i) + ". Error description: " +
              re.what() + " Please check also the log messages for details.";
      throw;
    }
    difc[i] = alg->getProperty("Difc");
    tzero[i] = alg->getProperty("Zero");

    g_log.notice() << "Bank " << i + 1 << " calibrated, "
                   << "difc: " << difc[i] << ", zero: " << tzero[i]
                   << std::endl;
  }

  // Double horror: 1st use a python script
  // 2nd: because runPythonCode does this by emitting a signal that goes to
  // MantidPlot,
  // it has to be done in the view (which is a UserSubWindow).
  Poco::Path outFullPath(cs.m_inputDirCalib);
  outFullPath.append(outFilename);
  m_view->writeOutCalibFile(outFullPath.toString(), difc, tzero);
  g_log.notice() << "Calibration file written as " << outFullPath.toString()
                 << std::endl;
}

/**
 * Parses the name of a calibration file and guesses the instrument,
 * vanadium and ceria run numbers, assuming that the name has been
 * build with buildCalibrateSuggestedFilename().
 *
 * @todo this is currently based on the filename. This method should
 * do a basic sanity check that the file is actually a calibration
 * file (IPARM file for GSAS), and get the numbers from the file not
 * from the name of the file. Leaving this as a TODO for now, as the
 * file template and contents are still subject to change.
 *
 * @param path full path to a calibration file
 * @param instName instrument used in the file
 * @param vanNo number of the Vanadium run used for this calibration file
 * @param ceriaNo number of the Ceria run
 *
 * @throws invalid_argument with an informative message if tha name does
 * not look correct (does not follow the conventions).
 */
void EnggDiffractionPresenter::parseCalibrateFilename(const std::string &path,
                                                      std::string &instName,
                                                      std::string &vanNo,
                                                      std::string &ceriaNo) {
  instName = "";
  vanNo = "";
  ceriaNo = "";

  Poco::Path fullPath(path);
  const std::string filename = fullPath.getFileName();
  if (filename.empty()) {
    return;
  }

  const std::string explMsg =
      "Expected a file name like 'INSTR_vanNo_ceriaNo_....par', "
      "where INSTR is the instrument name and vanNo and ceriaNo are the "
      "numbers of the Vanadium and calibration sample (Ceria, CeO2) runs.";
  std::vector<std::string> parts;
  boost::split(parts, filename, boost::is_any_of("_"));
  if (parts.size() < 4) {
    throw std::invalid_argument(
        "Failed to find at least the 4 required parts of the file name.\n\n" +
        explMsg);
  }

  // check the rules on the file name
  if (g_enginxStr != parts[0]) {
    throw std::invalid_argument("The first component of the file name is not "
                                "the expected instrument name: " +
                                g_enginxStr + ".\n\n" + explMsg);
  }
  const std::string castMsg =
      "It is not possible to interpret as an integer number ";
  try {
    boost::lexical_cast<int>(parts[1]);
  } catch (std::runtime_error &) {
    throw std::invalid_argument(
        castMsg + "the Vanadium number part of the file name.\n\n" + explMsg);
  }
  try {
    boost::lexical_cast<int>(parts[2]);
  } catch (std::runtime_error &) {
    throw std::invalid_argument(
        castMsg + "the Ceria number part of the file name.\n\n" + explMsg);
  }

  instName = parts[0];
  vanNo = parts[1];
  ceriaNo = parts[2];
}

/**
 * Build a suggested name for a new calibration, by appending instrument name,
 * relevant run numbers, etc., like: ENGINX_241391_236516_both_banks.par
 *
 * @param vanNo number of the Vanadium run
 * @param ceriaNo number of the Ceria run
 *
 * @return Suggested name for a new calibration file, following
 * ENGIN-X practices
 */
std::string EnggDiffractionPresenter::buildCalibrateSuggestedFilename(
    const std::string &vanNo, const std::string &ceriaNo) const {
  // default and only one supported
  std::string instStr = g_enginxStr;
  std::string nameAppendix = "_both_banks";
  if ("ENGIN-X" != m_view->currentInstrument()) {
    instStr = "UNKNOWNINST";
    nameAppendix = "_calibration";
  }

  // default extension for calibration files
  const std::string calibExt = ".prm";
  std::string sugg =
      instStr + "_" + vanNo + "_" + ceriaNo + nameAppendix + calibExt;

  return sugg;
}

/**
 * Load precalculated results from Vanadium corrections previously
 * calculated.
 *
 * @param preIntegFilename filename (can be full path) where the
 * vanadium spectra integration table should be loaded from
 *
 * @param preCurvesFilename filename (can be full path) where the
 * vanadium per-bank curves should be loaded from
 *
 * @param vanIntegWS output (matrix) workspace loaded from the
 * precalculated Vanadium correction file, with the integration
 * resutls
 *
 * @param vanCurvesWS output (matrix) workspace loaded from the
 * precalculated Vanadium correction file, with the per-bank curves
 */
void EnggDiffractionPresenter::loadVanadiumPrecalcWorkspaces(
    const std::string &preIntegFilename, const std::string &preCurvesFilename,
    ITableWorkspace_sptr &vanIntegWS, MatrixWorkspace_sptr &vanCurvesWS) {
  AnalysisDataServiceImpl &ADS = Mantid::API::AnalysisDataService::Instance();

  auto alg = Algorithm::fromString("LoadNexus");
  alg->initialize();
  alg->setPropertyValue("Filename", preIntegFilename);
  std::string integWSName = "engggui_vanadium_integration_ws";
  alg->setPropertyValue("OutputWorkspace", integWSName);
  alg->execute();
  // alg->getProperty("OutputWorkspace");
  vanIntegWS = ADS.retrieveWS<ITableWorkspace>(integWSName);

  auto algCurves = Algorithm::fromString("LoadNexus");
  algCurves->initialize();
  algCurves->setPropertyValue("Filename", preCurvesFilename);
  std::string curvesWSName = "engggui_vanadium_curves_ws";
  algCurves->setPropertyValue("OutputWorkspace", curvesWSName);
  algCurves->execute();
  // algCurves->getProperty("OutputWorkspace");
  vanCurvesWS = ADS.retrieveWS<MatrixWorkspace>(curvesWSName);
}

void EnggDiffractionPresenter::calcVanadiumWorkspaces(
    const std::string &vanNo, ITableWorkspace_sptr &vanIntegWS,
    MatrixWorkspace_sptr &vanCurvesWS) {

  auto load = Algorithm::fromString("Load");
  load->initialize();
  load->setPropertyValue("Filename",
                         vanNo); // TODO more specific build Vanadium filename
  std::string vanWSName = "engggui_vanadium_ws";
  load->setPropertyValue("OutputWorkspace", vanWSName);
  load->execute();
  AnalysisDataServiceImpl &ADS = Mantid::API::AnalysisDataService::Instance();
  MatrixWorkspace_sptr vanWS = ADS.retrieveWS<MatrixWorkspace>(vanWSName);
  // TODO: maybe use setChild() and then load->getProperty("OutputWorkspace");

  auto alg = Algorithm::fromString("EnggVanadiumCorrections");
  alg->initialize();
  alg->setProperty("VanadiumWorkspace", vanWS);
  std::string integName = "engggui_van_integration_ws";
  alg->setPropertyValue("IntegrationWorkspace", integName);
  std::string curvesName = "engggui_van_curves_ws";
  alg->setPropertyValue("CurvesWorkspace", curvesName);
  alg->execute();

  ADS.remove(vanWSName);

  vanIntegWS = ADS.retrieveWS<ITableWorkspace>(integName);
  vanCurvesWS = ADS.retrieveWS<MatrixWorkspace>(curvesName);
}

/**
 * builds the expected names of the precalculated Vanadium correction
 * files and tells if both files are found, similar to:
 * ENGINX_precalculated_vanadium_run000236516_integration.nxs
 * ENGINX_precalculated_vanadium_run00236516_bank_curves.nxs
 *
 * @param vanNo
 * @param inputDirCalib
 * @param preIntegFilename if not found on disk, the string is set as empty
 * @param preCurvesFilename if not found on disk, the string is set as empty
 * @param found true if both files are found and (re-)usable
 */
void EnggDiffractionPresenter::findPrecalcVanadiumCorrFilenames(
    const std::string &vanNo, const std::string &inputDirCalib,
    std::string &preIntegFilename, std::string &preCurvesFilename,
    bool &found) {
  found = false;

  const std::string runNo = std::string(2, '0').append(vanNo);
  preIntegFilename =
      g_enginxStr + "_precalculated_vanadium_run" + runNo + "_integration.nxs";

  preCurvesFilename =
      g_enginxStr + "_precalculated_vanadium_run" + runNo + "_bank_curves.nxs";

  Poco::Path pathInteg(inputDirCalib);
  pathInteg.append(preIntegFilename);

  Poco::Path pathCurves(inputDirCalib);
  pathCurves.append(preCurvesFilename);

  if (Poco::File(pathInteg).exists() && Poco::File(pathCurves).exists()) {
    preIntegFilename = pathInteg.toString();
    preCurvesFilename = pathCurves.toString();
    found = true;
  }
}

/**
 * Produce the two workspaces that are required to apply Vanadium
 * corrections. Try to load them if precalculated results are
 * available from files, otherwise load the source Vanadium run
 * workspace and do the calculations.
 *
 * @param vanNo Vanadium run number
 *
 * @param inputDirCalib The 'calibration files' input directory given
 * in settings
 *
 * @param vanIntegWS workspace where to create/load the Vanadium
 * spectra integration
 *
 * @param vanCurvesWS workspace where to create/load the Vanadium
 * aggregated per-bank curve
 *
 * @param forceRecalc whether to calculate Vanadium corrections even
 * if the files of pre-calculated results are found
 */
void EnggDiffractionPresenter::loadOrCalcVanadiumWorkspaces(
    const std::string &vanNo, const std::string &inputDirCalib,
    ITableWorkspace_sptr &vanIntegWS, MatrixWorkspace_sptr &vanCurvesWS,
    bool forceRecalc) {
  bool foundPrecalc = false;

  std::string preIntegFilename, preCurvesFilename;
  findPrecalcVanadiumCorrFilenames(vanNo, inputDirCalib, preIntegFilename,
                                   preCurvesFilename, foundPrecalc);

  if (forceRecalc || !foundPrecalc) {
    g_log.notice()
        << "Calculating Vanadium corrections. This may take a few seconds..."
        << std::endl;
    try {
      calcVanadiumWorkspaces(vanNo, vanIntegWS, vanCurvesWS);
    } catch (std::invalid_argument &ia) {
      g_log.error() << "Failed to calculate Vanadium corrections. "
                       "There was an error in the execution of the algorithms "
                       "required to calculate Vanadium corrections. Some "
                       "properties passed to the algorithms were invalid. "
                       "This is possibly because some of the settings are not "
                       "consistent. Please check the log messages for "
                       "details. Details: " +
                           std::string(ia.what()) << std::endl;
      throw;
    } catch (std::runtime_error &re) {
      g_log.error() << "Failed to calculate Vanadium corrections. "
                       "There was an error while executing one of the "
                       "algorithms used to perform Vanadium corrections. "
                       "There was no obvious error in the input properties "
                       "but the algorithm failed. Please check the log "
                       "messages for details." +
                           std::string(re.what()) << std::endl;
      throw;
    }
  } else {
    g_log.notice()
        << "Found precalculated Vanadium correction features for Vanadium run "
        << vanNo << ". Re-using these files: " << preIntegFilename << ", and "
        << preCurvesFilename << std::endl;
    try {
      loadVanadiumPrecalcWorkspaces(preIntegFilename, preCurvesFilename,
                                    vanIntegWS, vanCurvesWS);
    } catch (std::invalid_argument &ia) {
      g_log.error() << "Error while loading precalculated Vanadium corrections",
          "The files with precalculated Vanadium corection features (spectra "
          "integration and per-bank curves) were found (with names '" +
              preIntegFilename + "' and '" + preCurvesFilename +
              "', respectively, but there was a problem with the inputs to the "
              "load algorithms to load them: " +
              std::string(ia.what());
      throw;
    } catch (std::runtime_error &re) {
      g_log.error() << "Error while loading precalculated Vanadium corrections",
          "The files with precalculated Vanadium corection features (spectra "
          "integration and per-bank curves) were found (with names '" +
              preIntegFilename + "' and '" + preCurvesFilename +
              "', respectively, but there was a problem while loading them. "
              "Please check the log messages for details. You might want to "
              "delete those files or force recalculations (in settings). Error "
              "details: " +
              std::string(re.what());
      throw;
    }
  }
}

} // namespace CustomInterfaces
} // namespace MantidQt
