#ifndef MANTID_ALGORITHMS_FILTEREVENTSTEST_H_
#define MANTID_ALGORITHMS_FILTEREVENTSTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidAPI/DetectorInfo.h"
#include "MantidAPI/SpectrumInfo.h"
#include "MantidAPI/TableRow.h"
#include "MantidAlgorithms/FilterEvents.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidDataObjects/EventList.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidDataObjects/Events.h"
#include "MantidDataObjects/SplittersWorkspace.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidKernel/PhysicalConstants.h"
#include "MantidKernel/TimeSeriesProperty.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"

#include <random>

using namespace Mantid;
using namespace Mantid::Algorithms;
using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::Geometry;

using namespace std;

class FilterEventsTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static FilterEventsTest *createSuite() { return new FilterEventsTest(); }
  static void destroySuite(FilterEventsTest *suite) { delete suite; }

  //----------------------------------------------------------------------------------------------
  /** Test initialization
    */
  void test_Initialization() {
    FilterEvents alg;
    alg.initialize();

    TS_ASSERT(alg.isInitialized());
  }

  //----------------------------------------------------------------------------------------------
  /** Test create event workspace and splitters
    * In all the tests below:
    * (1) 10 detectors
    * (2) Run starts @ 20000000000 seconds
    * (3) Pulse length = 100*1000*1000 seconds
    * (4) Within one pulse, two consecutive events/neutrons is apart for
   * 10*1000*1000 seconds
    * (5) "Experiment": 5 pulse times.  10 events in each pulse
    */
  void test_CreatedEventWorskpaceAndSplitter() {
    int64_t runstart_i64 = 20000000000;
    int64_t pulsedt = 100 * 1000 * 1000;
    int64_t tofdt = 10 * 1000 * 1000;
    size_t numpulses = 5;

    EventWorkspace_sptr eventws =
        createEventWorkspace(runstart_i64, pulsedt, tofdt, numpulses);

    TS_ASSERT_EQUALS(eventws->getNumberEvents(), 500);

    EventList elist = eventws->getSpectrum(0);
    TS_ASSERT_EQUALS(elist.getNumberEvents(), 50);
    TS_ASSERT(elist.hasDetectorID(1));

    SplittersWorkspace_sptr splittersws =
        createSplitter(runstart_i64, pulsedt, tofdt);
    TS_ASSERT_EQUALS(splittersws->getNumberSplitters(), 5);

    return;
  }

  //----------------------------------------------------------------------------------------------
  /**  Filter events without any correction
    *  Event workspace:
    * (1) 10 detectors
    * (2) Run starts @ 20000000000 seconds
    * (3) Pulse length = 100*1000*1000 seconds
    * (4) Within one pulse, two consecutive events/neutrons is apart for
   *10*1000*1000 seconds
    * (5) "Experiment": 5 pulse times.  10 events in each pulse
    *
    * In this test
   *  (1) Leave correction table workspace empty
   *  (2) Count events in each output including "-1", the excluded/unselected
   *events
   */
  void test_FilterNoCorrection() {
    // Create EventWorkspace and SplittersWorkspace
    int64_t runstart_i64 = 20000000000;
    int64_t pulsedt = 100 * 1000 * 1000;
    int64_t tofdt = 10 * 1000 * 1000;
    size_t numpulses = 5;

    EventWorkspace_sptr inpWS =
        createEventWorkspace(runstart_i64, pulsedt, tofdt, numpulses);
    AnalysisDataService::Instance().addOrReplace("Test02", inpWS);

    SplittersWorkspace_sptr splws =
        createSplitter(runstart_i64, pulsedt, tofdt);
    AnalysisDataService::Instance().addOrReplace("Splitter02", splws);

    FilterEvents filter;
    filter.initialize();

    // Set properties
    filter.setProperty("InputWorkspace", "Test02");
    filter.setProperty("OutputWorkspaceBaseName", "FilteredWS01");
    filter.setProperty("SplitterWorkspace", "Splitter02");
    filter.setProperty("OutputTOFCorrectionWorkspace", "CorrectionWS");

    // Execute
    TS_ASSERT_THROWS_NOTHING(filter.execute());
    TS_ASSERT(filter.isExecuted());

    // Get output
    int numsplittedws = filter.getProperty("NumberOutputWS");
    TS_ASSERT_EQUALS(numsplittedws, 4);

    // Check Workspace group 0
    EventWorkspace_sptr filteredws0 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS01_0"));
    TS_ASSERT(filteredws0);
    TS_ASSERT_EQUALS(filteredws0->getNumberHistograms(), 10);
    TS_ASSERT_EQUALS(filteredws0->getSpectrum(0).getNumberEvents(), 4);
    TS_ASSERT_EQUALS(filteredws0->run().getProtonCharge(), 10);

    // Check Workspace group 1
    EventWorkspace_sptr filteredws1 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS01_1"));
    TS_ASSERT(filteredws1);
    TS_ASSERT_EQUALS(filteredws1->getSpectrum(1).getNumberEvents(), 16);
    TS_ASSERT_EQUALS(filteredws1->run().getProtonCharge(), 11);

    // Check Workspace group 2
    EventWorkspace_sptr filteredws2 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS01_2"));
    TS_ASSERT(filteredws2);
    TS_ASSERT_EQUALS(filteredws2->getSpectrum(1).getNumberEvents(), 21);
    TS_ASSERT_EQUALS(filteredws2->run().getProtonCharge(), 21);

    EventList elist3 = filteredws2->getSpectrum(3);
    elist3.sortPulseTimeTOF();

    TofEvent eventmin = elist3.getEvent(0);
    TS_ASSERT_EQUALS(eventmin.pulseTime().totalNanoseconds(),
                     runstart_i64 + pulsedt * 2);
    TS_ASSERT_DELTA(eventmin.tof(), 0, 1.0E-4);

    TofEvent eventmax = elist3.getEvent(20);
    TS_ASSERT_EQUALS(eventmax.pulseTime().totalNanoseconds(),
                     runstart_i64 + pulsedt * 4);
    TS_ASSERT_DELTA(eventmax.tof(), static_cast<double>(tofdt * 6 / 1000),
                    1.0E-4);

    // Clean up
    AnalysisDataService::Instance().remove("Test02");
    AnalysisDataService::Instance().remove("Splitter02");
    std::vector<std::string> outputwsnames =
        filter.getProperty("OutputWorkspaceNames");
    for (size_t i = 0; i < outputwsnames.size(); ++i) {
      AnalysisDataService::Instance().remove(outputwsnames[i]);
    }

    return;
  }

  //----------------------------------------------------------------------------------------------
  /**  Filter events without any correction and test for user-specified
   *workspace starting value
    *  Event workspace:
    * (1) 10 detectors
    * (2) Run starts @ 20000000000 seconds
    * (3) Pulse length = 100*1000*1000 seconds
    * (4) Within one pulse, two consecutive events/neutrons is apart for
   *10*1000*1000 seconds
    * (5) "Experiment": 5 pulse times.  10 events in each pulse
    *
    * In this test
   *  (1) Leave correction table workspace empty
   *  (2) Count events in each output including "-1", the excluded/unselected
   *events
   */
  void Passed_test_FilterWOCorrection2() {
    // Create EventWorkspace and SplittersWorkspace
    int64_t runstart_i64 = 20000000000;
    int64_t pulsedt = 100 * 1000 * 1000;
    int64_t tofdt = 10 * 1000 * 1000;
    size_t numpulses = 5;

    EventWorkspace_sptr inpWS =
        createEventWorkspace(runstart_i64, pulsedt, tofdt, numpulses);
    AnalysisDataService::Instance().addOrReplace("Test02", inpWS);

    SplittersWorkspace_sptr splws =
        createSplitter(runstart_i64, pulsedt, tofdt);
    AnalysisDataService::Instance().addOrReplace("Splitter02", splws);

    FilterEvents filter;
    filter.initialize();

    // Set properties
    filter.setProperty("InputWorkspace", "Test02");
    filter.setProperty("OutputWorkspaceBaseName", "FilteredWS01");
    filter.setProperty("SplitterWorkspace", "Splitter02");
    filter.setProperty("OutputWorkspaceIndexedFrom1", true);

    // Execute
    TS_ASSERT_THROWS_NOTHING(filter.execute());
    TS_ASSERT(filter.isExecuted());

    // Get output
    int numsplittedws = filter.getProperty("NumberOutputWS");
    TS_ASSERT_EQUALS(numsplittedws, 3);

    // 4.1 Workspace group 0
    EventWorkspace_sptr filteredws0 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS01_1"));
    TS_ASSERT(filteredws0);
    TS_ASSERT_EQUALS(filteredws0->getNumberHistograms(), 10);
    TS_ASSERT_EQUALS(filteredws0->getSpectrum(0).getNumberEvents(), 4);

    // 4.2 Workspace group 1
    EventWorkspace_sptr filteredws1 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS01_2"));
    TS_ASSERT(filteredws1);
    TS_ASSERT_EQUALS(filteredws1->getSpectrum(1).getNumberEvents(), 16);

    // 4.3 Workspace group 2
    EventWorkspace_sptr filteredws2 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS01_3"));
    TS_ASSERT(filteredws2);
    TS_ASSERT_EQUALS(filteredws2->getSpectrum(1).getNumberEvents(), 21);

    EventList elist3 = filteredws2->getSpectrum(3);
    elist3.sortPulseTimeTOF();

    TofEvent eventmin = elist3.getEvent(0);
    TS_ASSERT_EQUALS(eventmin.pulseTime().totalNanoseconds(),
                     runstart_i64 + pulsedt * 2);
    TS_ASSERT_DELTA(eventmin.tof(), 0, 1.0E-4);

    TofEvent eventmax = elist3.getEvent(20);
    TS_ASSERT_EQUALS(eventmax.pulseTime().totalNanoseconds(),
                     runstart_i64 + pulsedt * 4);
    TS_ASSERT_DELTA(eventmax.tof(), static_cast<double>(tofdt * 6 / 1000),
                    1.0E-4);

    // 5. Clean up
    AnalysisDataService::Instance().remove("Test02");
    AnalysisDataService::Instance().remove("Splitter02");
    std::vector<std::string> outputwsnames =
        filter.getProperty("OutputWorkspaceNames");
    for (size_t i = 0; i < outputwsnames.size(); ++i) {
      AnalysisDataService::Instance().remove(outputwsnames[i]);
    }

    return;
  }

  //----------------------------------------------------------------------------------------------
  /**  Filter test with TOF correction
    */
  void Passed_test_FilterWithCustumizedCorrection() {
    // 1. Create EventWorkspace and SplittersWorkspace
    int64_t runstart_i64 = 20000000000;
    int64_t pulsedt = 100 * 1000 * 1000;
    int64_t tofdt = 10 * 1000 * 1000;
    size_t numpulses = 5;

    EventWorkspace_sptr inpWS =
        createEventWorkspace(runstart_i64, pulsedt, tofdt, numpulses);
    AnalysisDataService::Instance().addOrReplace("EventData", inpWS);

    SplittersWorkspace_sptr splws =
        createFastFreqLogSplitter(runstart_i64, pulsedt, tofdt, numpulses);
    AnalysisDataService::Instance().addOrReplace("SplitterTableX", splws);
    TS_ASSERT_EQUALS(splws->rowCount(), static_cast<size_t>(numpulses) * 2);

    TableWorkspace_sptr timecorrws = createTimeCorrectionTable(inpWS);
    AnalysisDataService::Instance().addOrReplace("TimeCorrectionTableX",
                                                 timecorrws);
    TS_ASSERT_EQUALS(timecorrws->rowCount(), inpWS->getNumberHistograms());

    FilterEvents filter;
    filter.initialize();

    // 2. Set properties
    TS_ASSERT_THROWS_NOTHING(filter.setProperty("InputWorkspace", "EventData"));
    TS_ASSERT_THROWS_NOTHING(
        filter.setProperty("OutputWorkspaceBaseName", "SplittedDataX"));
    TS_ASSERT_THROWS_NOTHING(
        filter.setProperty("CorrectionToSample", "Customized"));
    TS_ASSERT_THROWS_NOTHING(filter.setProperty(
        "DetectorTOFCorrectionWorkspace", "TimeCorrectionTableX"));
    TS_ASSERT_THROWS_NOTHING(filter.setProperty("SplitterWorkspace", splws));

    // 3. Execute
    TS_ASSERT_THROWS_NOTHING(filter.execute());
    TS_ASSERT(filter.isExecuted());

    // 4. Get output
    // 4.1 Workspace group 0
    EventWorkspace_sptr filteredws0 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("SplittedDataX_0"));
    TS_ASSERT(filteredws0);
    TS_ASSERT_EQUALS(filteredws0->getNumberHistograms(), 10);
    TS_ASSERT_EQUALS(filteredws0->getSpectrum(0).getNumberEvents(), 15);
    TS_ASSERT_EQUALS(filteredws0->getSpectrum(9).getNumberEvents(), 15);
    TS_ASSERT_EQUALS(filteredws0->run().getProtonCharge(), 5);

    // 4.2 Workspace group 1
    EventWorkspace_sptr filteredws1 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("SplittedDataX_1"));
    TS_ASSERT(filteredws1);
    TS_ASSERT_EQUALS(filteredws1->getSpectrum(1).getNumberEvents(), 10);
    TS_ASSERT_EQUALS(filteredws0->run().getProtonCharge(), 5);

    // 4.3 Some individual events
    EventList elist3 = filteredws1->getSpectrum(3);
    elist3.sortPulseTimeTOF();

    if (elist3.getNumberEvents() > 0) {
      TofEvent eventmin = elist3.getEvent(0);
      TS_ASSERT_EQUALS(eventmin.pulseTime().totalNanoseconds(), runstart_i64);
      TS_ASSERT_DELTA(eventmin.tof(), 80 * 1000, 1.0E-4);
    }

    // 5. Clean
    AnalysisDataService::Instance().remove("EventData");
    AnalysisDataService::Instance().remove("TimeCorrectionTableX");
    AnalysisDataService::Instance().remove("SplitterTableX");

    std::vector<std::string> outputwsnames =
        filter.getProperty("OutputWorkspaceNames");
    for (size_t i = 0; i < outputwsnames.size(); ++i) {
      AnalysisDataService::Instance().remove(outputwsnames[i]);
    }

    return;
  }

  //----------------------------------------------------------------------------------------------
  /** Test filtering with correction of direct geometry
    */
  void Passed_test_FilterElasticCorrection() {
    EventWorkspace_sptr ws = createEventWorkspaceElastic(0, 1000000);
    AnalysisDataService::Instance().addOrReplace("MockElasticEventWS", ws);
    TS_ASSERT_EQUALS(ws->getNumberEvents(), 10000);

    MatrixWorkspace_sptr splws = createMatrixSplittersElastic();
    AnalysisDataService::Instance().addOrReplace("SplitterTableX", splws);

    // Run the filtering
    FilterEvents filter;
    filter.initialize();

    filter.setPropertyValue("InputWorkspace", "MockElasticEventWS");
    filter.setProperty("OutputWorkspaceBaseName", "SplittedDataElastic");
    filter.setProperty("CorrectionToSample", "Elastic");
    filter.setProperty("SplitterWorkspace", "SplitterTableX");

    TS_ASSERT_THROWS_NOTHING(filter.execute());
    TS_ASSERT(filter.isExecuted());

    // Check number of output workspaces
    std::vector<std::string> vecwsname =
        filter.getProperty("OutputWorkspaceNames");
    TS_ASSERT_EQUALS(vecwsname.size(), 9);

    EventWorkspace_sptr ws5 = boost::dynamic_pointer_cast<EventWorkspace>(
        AnalysisDataService::Instance().retrieve("SplittedDataElastic_5"));
    TS_ASSERT(ws5);
    if (ws5) {
      TS_ASSERT_EQUALS(ws5->getNumberEvents(), 0);
    }

    EventWorkspace_sptr ws7 = boost::dynamic_pointer_cast<EventWorkspace>(
        AnalysisDataService::Instance().retrieve("SplittedDataElastic_7"));
    TS_ASSERT(ws7);
    if (ws7) {
      TS_ASSERT_EQUALS(ws7->getNumberEvents(), 10);
    }

    // Check individual events
    EventList &ev0 = ws7->getSpectrum(0);
    TS_ASSERT_EQUALS(ev0.getNumberEvents(), 1);
    std::vector<double> vectofs = ev0.getTofs();
    TS_ASSERT_DELTA(vectofs[0], 272.0, 0.001);

    // Delete all the workspaces generated here
    AnalysisDataService::Instance().remove("MockDirectEventWS");
    AnalysisDataService::Instance().remove("SplitterTableX");
    for (size_t i = 0; i < vecwsname.size(); ++i) {
      AnalysisDataService::Instance().remove(vecwsname[i]);
    }

    return;
  }

  //----------------------------------------------------------------------------------------------
  /** Test filtering with correction of direct geometry
    */
  void Passed_test_FilterDGCorrection() {
    EventWorkspace_sptr ws = createEventWorkspaceDirect(0, 1000000);
    AnalysisDataService::Instance().addOrReplace("MockDirectEventWS", ws);

    MatrixWorkspace_sptr splws = createMatrixSplittersDG();
    AnalysisDataService::Instance().addOrReplace("SplitterTableX", splws);

    // Run the filtering
    FilterEvents filter;
    filter.initialize();

    filter.setProperty("InputWorkspace", ws->getName());
    filter.setProperty("OutputWorkspaceBaseName", "SplittedDataDG");
    filter.setProperty("CorrectionToSample", "Direct");
    filter.setProperty("SplitterWorkspace", "SplitterTableX");

    TS_ASSERT_THROWS_NOTHING(filter.execute());
    TS_ASSERT(filter.isExecuted());

    // Check
    EventWorkspace_sptr ws5 = boost::dynamic_pointer_cast<EventWorkspace>(
        AnalysisDataService::Instance().retrieve("SplittedDataDG_5"));
    TS_ASSERT(ws5);
    if (ws5) {
      TS_ASSERT_EQUALS(ws5->getNumberEvents(), 0);
    }

    EventWorkspace_sptr ws7 = boost::dynamic_pointer_cast<EventWorkspace>(
        AnalysisDataService::Instance().retrieve("SplittedDataDG_7"));
    TS_ASSERT(ws7);
    if (ws7) {
      TS_ASSERT_EQUALS(ws7->getNumberEvents(), ws7->getNumberHistograms());
    }

    // FIXME - Should find a way to delete all workspaces holding splitted
    // events

    AnalysisDataService::Instance().remove("MockDirectEventWS");
    AnalysisDataService::Instance().remove("SplitterTableX");
    std::vector<std::string> outputwsnames =
        filter.getProperty("OutputWorkspaceNames");
    for (size_t i = 0; i < outputwsnames.size(); ++i) {
      AnalysisDataService::Instance().remove(outputwsnames[i]);
    }

    return;
  }

  //----------------------------------------------------------------------------------------------
  /** Test filtering with correction to indirect geometry inelastic instrument
    */
  void Passed_test_FilterIndirectGeometryCorrection() {
    // Create workspaces for filtering
    EventWorkspace_sptr ws = createEventWorkspaceInDirect(0, 1000000);
    AnalysisDataService::Instance().addOrReplace("MockIndirectEventWS", ws);

    MatrixWorkspace_sptr splws = createMatrixSplittersDG();
    AnalysisDataService::Instance().addOrReplace("SplitterTableX", splws);

    // Run the filtering
    FilterEvents filter;
    filter.initialize();

    filter.setProperty("InputWorkspace", "MockIndirectEventWS");
    filter.setProperty("OutputWorkspaceBaseName", "SplittedDataDG");
    filter.setProperty("CorrectionToSample", "Indirect");
    filter.setProperty("SplitterWorkspace", "SplitterTableX");
    filter.setProperty("OutputTOFCorrectionWorkspace", "MockIndGeoCorrWS");

    TS_ASSERT_THROWS_NOTHING(filter.execute());
    TS_ASSERT(filter.isExecuted());

    // Check
    MatrixWorkspace_sptr outcorrws =
        boost::dynamic_pointer_cast<MatrixWorkspace>(
            AnalysisDataService::Instance().retrieve("MockIndGeoCorrWS"));
    TS_ASSERT(outcorrws);
    if (outcorrws) {
      TS_ASSERT_EQUALS(outcorrws->getNumberHistograms(),
                       ws->getNumberHistograms());
      TS_ASSERT_EQUALS(outcorrws->x(0).size(), 2);

      const auto &spectrumInfo = ws->spectrumInfo();

      for (size_t iws = 0; iws < outcorrws->getNumberHistograms(); ++iws) {
        const ParameterMap &pmap = ws->constInstrumentParameters();

        const auto &det = spectrumInfo.detector(iws);
        Parameter_sptr par = pmap.getRecursive(&det, "Efixed");
        double efix = par->value<double>();

        double l2 = spectrumInfo.l2(iws);

        double shift = -l2 / sqrt(efix * 2. * PhysicalConstants::meV /
                                  PhysicalConstants::NeutronMass);

        TS_ASSERT_DELTA(outcorrws->y(iws)[0], 1., 1.0E-9);
        TS_ASSERT_DELTA(outcorrws->y(iws)[1], shift, 1.0E-9);
      }
    }

    // Clean
    AnalysisDataService::Instance().remove("MockIndirectEventWS");
    std::vector<std::string> outputwsnames =
        filter.getProperty("OutputWorkspaceNames");
    for (size_t i = 0; i < outputwsnames.size(); ++i) {
      AnalysisDataService::Instance().remove(outputwsnames[i]);
    }

    return;
  }

  //----------------------------------------------------------------------------------------------
  /**  Filter events without any correction and test for splitters in
   *MatrixWorkspace format
   *   and the time given for splitters is relative
   *  Event workspace:
   * (1) 10 detectors
    * (2) Run starts @ 20000000000 seconds
    * (3) Pulse length = 100*1000*1000 seconds
    * (4) Within one pulse, two consecutive events/neutrons is apart for
   *10*1000*1000 seconds
    * (5) "Experiment": 5 pulse times.  10 events in each pulse
    *
    * In this test
   *  (1) Leave correction table workspace empty
   *  (2) Count events in each output including "-1", the excluded/unselected
   *events
   */
  void test_FilterRelativeTime() {
    // Create EventWorkspace and SplittersWorkspace
    int64_t runstart_i64 = 20000000000;
    int64_t pulsedt = 100 * 1000 * 1000;
    int64_t tofdt = 10 * 1000 * 1000;
    size_t numpulses = 5;

    EventWorkspace_sptr inpWS =
        createEventWorkspace(runstart_i64, pulsedt, tofdt, numpulses);
    AnalysisDataService::Instance().addOrReplace("Test10", inpWS);

    API::MatrixWorkspace_sptr splws = createMatrixSplitter(0, pulsedt, tofdt);
    AnalysisDataService::Instance().addOrReplace("Splitter10", splws);

    FilterEvents filter;
    filter.initialize();

    // Set properties
    filter.setProperty("InputWorkspace", "Test10");
    filter.setProperty("OutputWorkspaceBaseName", "FilteredWS10");
    filter.setProperty("SplitterWorkspace", "Splitter10");
    filter.setProperty("RelativeTime", true);
    filter.setProperty("OutputWorkspaceIndexedFrom1", false);

    // Execute
    TS_ASSERT_THROWS_NOTHING(filter.execute());
    TS_ASSERT(filter.isExecuted());

    // Get 3 output workspaces
    int numsplittedws = filter.getProperty("NumberOutputWS");
    TS_ASSERT_EQUALS(numsplittedws, 3);

    std::vector<std::string> output_ws_vector =
        filter.getProperty("OutputWorkspaceNames");
    for (size_t i = 0; i < output_ws_vector.size(); ++i)
      std::cout << "Output workspace " << i << ": " << output_ws_vector[i]
                << "\n";

    // Workspace 0
    EventWorkspace_sptr filteredws0 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS10_0"));
    TS_ASSERT(filteredws0);
    TS_ASSERT_EQUALS(filteredws0->getNumberHistograms(), 10);
    TS_ASSERT_EQUALS(filteredws0->getSpectrum(0).getNumberEvents(), 3);

    // Workspace 1
    EventWorkspace_sptr filteredws1 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS10_1"));
    TS_ASSERT(filteredws1);
    TS_ASSERT_EQUALS(filteredws1->getSpectrum(1).getNumberEvents(), 16);

    // Workspace 2
    EventWorkspace_sptr filteredws2 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS10_2"));
    TS_ASSERT(filteredws2);
    TS_ASSERT_EQUALS(filteredws2->getSpectrum(1).getNumberEvents(), 27);

    // Check spectrum 3 of workspace 2
    EventList elist3 = filteredws2->getSpectrum(3);
    elist3.sortPulseTimeTOF();

    TofEvent eventmin = elist3.getEvent(0);
    TS_ASSERT_EQUALS(eventmin.pulseTime().totalNanoseconds(),
                     runstart_i64 + pulsedt * 2);
    TS_ASSERT_DELTA(eventmin.tof(), 0, 1.0E-4);

    TofEvent eventmax = elist3.getEvent(26);
    TS_ASSERT_EQUALS(eventmax.pulseTime().totalNanoseconds(),
                     runstart_i64 + pulsedt * 4);
    TS_ASSERT_DELTA(eventmax.tof(), static_cast<double>(tofdt * 6 / 1000),
                    1.0E-4);

    //  Test the sample logs
    std::vector<std::string> outputwsnames =
        filter.getProperty("OutputWorkspaceNames");
    for (size_t i = 0; i < outputwsnames.size(); ++i) {
      EventWorkspace_sptr filtered_ws =
          boost::dynamic_pointer_cast<DataObjects::EventWorkspace>(
              AnalysisDataService::Instance().retrieve(outputwsnames[i]));

      TS_ASSERT(filtered_ws->run().hasProperty("LogA"));
      TS_ASSERT(filtered_ws->run().hasProperty("LogB"));
      TS_ASSERT(filtered_ws->run().hasProperty("LogC"));

      Kernel::Property *logA = filtered_ws->run().getProperty("LogA");
      std::string valueA = logA->value();
      TS_ASSERT_EQUALS(valueA.compare("A"), 0);

      TS_ASSERT(filtered_ws->run().hasProperty("slow_int_log"));
      Kernel::TimeSeriesProperty<int> *intlog =
          dynamic_cast<Kernel::TimeSeriesProperty<int> *>(
              filtered_ws->run().getProperty("slow_int_log"));
      TS_ASSERT(intlog);
    }

    // clean up all the workspaces generated
    AnalysisDataService::Instance().remove("Test10");
    AnalysisDataService::Instance().remove("Splitter10");
    for (size_t i = 0; i < outputwsnames.size(); ++i) {
      AnalysisDataService::Instance().remove(outputwsnames[i]);
    }

    return;
  }

  //----------------------------------------------------------------------------------------------
  /**  Filter events without any correction and test for splitters in
   *    TableWorkspace filter format
   *    and the time given for splitters is relative
   *
   *  It is exacly the same as unit test: test_FilterRelativeTime()
   *
   *  Event workspace:
   * (1) 10 detectors
   * (2) Run starts @ 20000000000 seconds
   * (3) Pulse length = 100*1000*1000 seconds
   * (4) Within one pulse, two consecutive events/neutrons is apart for
   * 10*1000*1000 seconds
   * (5) "Experiment": 5 pulse times.  10 events in each pulse
   *
   * In this test
   *  (1) Leave correction table workspace empty
   *  (2) Count events in each output including "-1", the excluded/unselected
   * events
   */
  void Passed_test_tableSplitter() {
    // Create EventWorkspace and SplittersWorkspace
    int64_t runstart_i64 = 20000000000;
    int64_t pulsedt = 100 * 1000 * 1000;
    int64_t tofdt = 10 * 1000 * 1000;
    size_t numpulses = 5;

    EventWorkspace_sptr inpWS =
        createEventWorkspace(runstart_i64, pulsedt, tofdt, numpulses);
    AnalysisDataService::Instance().addOrReplace("Test11", inpWS);

    DataObjects::TableWorkspace_sptr splws =
        createTableSplitters(0, pulsedt, tofdt);
    AnalysisDataService::Instance().addOrReplace("TableSplitter1", splws);

    FilterEvents filter;
    filter.initialize();

    // Set properties
    filter.setProperty("InputWorkspace", "Test11");
    filter.setProperty("OutputWorkspaceBaseName", "FilteredWS_FromTable");
    filter.setProperty("SplitterWorkspace", "TableSplitter1");
    filter.setProperty("RelativeTime", true);
    filter.setProperty("OutputWorkspaceIndexedFrom1", true);

    // Execute
    TS_ASSERT_THROWS_NOTHING(filter.execute());
    TS_ASSERT(filter.isExecuted());

    // Get 3 output workspaces
    int numsplittedws = filter.getProperty("NumberOutputWS");
    TS_ASSERT_EQUALS(numsplittedws, 3);

    std::vector<std::string> output_ws_vector =
        filter.getProperty("OutputWorkspaceNames");
    for (size_t i = 0; i < output_ws_vector.size(); ++i)
      std::cout << "Output workspace " << i << ": " << output_ws_vector[i]
                << "\n";

    // Workspace 0
    EventWorkspace_sptr filteredws0 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS_FromTable_A"));
    TS_ASSERT(filteredws0);
    TS_ASSERT_EQUALS(filteredws0->getNumberHistograms(), 10);
    TS_ASSERT_EQUALS(filteredws0->getSpectrum(0).getNumberEvents(), 3);

    // Workspace 1
    EventWorkspace_sptr filteredws1 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS_FromTable_B"));
    TS_ASSERT(filteredws1);
    TS_ASSERT_EQUALS(filteredws1->getSpectrum(1).getNumberEvents(), 16);

    // Workspace 2
    EventWorkspace_sptr filteredws2 =
        boost::dynamic_pointer_cast<EventWorkspace>(
            AnalysisDataService::Instance().retrieve("FilteredWS_FromTable_C"));
    TS_ASSERT(filteredws2);
    TS_ASSERT_EQUALS(filteredws2->getSpectrum(1).getNumberEvents(), 27);

    // Check spectrum 3 of workspace 2
    EventList elist3 = filteredws2->getSpectrum(3);
    elist3.sortPulseTimeTOF();

    TofEvent eventmin = elist3.getEvent(0);
    TS_ASSERT_EQUALS(eventmin.pulseTime().totalNanoseconds(),
                     runstart_i64 + pulsedt * 2);
    TS_ASSERT_DELTA(eventmin.tof(), 0, 1.0E-4);

    TofEvent eventmax = elist3.getEvent(26);
    TS_ASSERT_EQUALS(eventmax.pulseTime().totalNanoseconds(),
                     runstart_i64 + pulsedt * 4);
    TS_ASSERT_DELTA(eventmax.tof(), static_cast<double>(tofdt * 6 / 1000),
                    1.0E-4);

    // 5. Clean up
    AnalysisDataService::Instance().remove("Test11");
    AnalysisDataService::Instance().remove("TableSplitter1");
    std::vector<std::string> outputwsnames =
        filter.getProperty("OutputWorkspaceNames");
    for (size_t i = 0; i < outputwsnames.size(); ++i) {
      AnalysisDataService::Instance().remove(outputwsnames[i]);
    }

    return;
  }

  //----------------------------------------------------------------------------------------------
  /** Create an EventWorkspace.  This workspace has
    * @param runstart_i64 : absolute run start time in int64_t format with unit
   * nanosecond
    * @param pulsedt : pulse length in int64_t format with unit nanosecond
    * @param todft : time interval between 2 adjacent event in same pulse in
   * int64_t format of unit nanosecond
    * @param numpulses : number of pulses in the event workspace
   */
  EventWorkspace_sptr createEventWorkspace(int64_t runstart_i64,
                                           int64_t pulsedt, int64_t tofdt,
                                           size_t numpulses) {
    // 1. Create an EventWorkspace with 10 detectors
    EventWorkspace_sptr eventWS =
        WorkspaceCreationHelper::createEventWorkspaceWithFullInstrument(10, 1,
                                                                        true);

    Kernel::DateAndTime runstart(runstart_i64);

    // 2. Set run_start time
    eventWS->mutableRun().addProperty("run_start", runstart.toISO8601String(),
                                      true);

    // create a pcharge log
    auto pchargeLog = Kernel::make_unique<Kernel::TimeSeriesProperty<double>>(
        "proton_charge");

    for (size_t i = 0; i < eventWS->getNumberHistograms(); i++) {
      auto &elist = eventWS->getSpectrum(i);

      for (int64_t pid = 0; pid < static_cast<int64_t>(numpulses); pid++) {
        int64_t pulsetime_i64 = pid * pulsedt + runstart.totalNanoseconds();
        Kernel::DateAndTime pulsetime(pulsetime_i64);
        pchargeLog->addValue(pulsetime, 1.);
        for (size_t e = 0; e < 10; e++) {
          double tof = static_cast<double>(e * tofdt / 1000);
          TofEvent event(tof, pulsetime);
          elist.addEventQuickly(event);
        }
      } // FOR each pulse
    }   // For each bank

    eventWS->mutableRun().addLogData(pchargeLog.release());
    eventWS->mutableRun().integrateProtonCharge();

    // add some arbitrary sample log for splitting or not splitting
    eventWS->mutableRun().addProperty(
        new Kernel::PropertyWithValue<std::string>("LogA", "A"));
    eventWS->mutableRun().addProperty(
        new Kernel::PropertyWithValue<std::string>("LogB", "B"));
    eventWS->mutableRun().addProperty(
        new Kernel::PropertyWithValue<std::string>("LogC", "C"), true);
    eventWS->mutableRun().addProperty(
        new Kernel::PropertyWithValue<std::string>("Title",
                                                   "Testing EventWorkspace"));

    // add an integer slow log
    auto int_tsp =
        Kernel::make_unique<Kernel::TimeSeriesProperty<int>>("slow_int_log");
    for (size_t i = 0; i < 10; ++i) {
      Kernel::DateAndTime log_time(runstart_i64 + 5 * pulsedt * i);
      int log_value = static_cast<int>(i + 1) * 20;
      int_tsp->addValue(log_time, log_value);
    }
    eventWS->mutableRun().addLogData(int_tsp.release());

    return eventWS;
  }

  //----------------------------------------------------------------------------------------------
  /** Create an EventWorkspace to mimic direct inelastic scattering insturment.
    * This workspace will have the same neutron events as the test case in
    *EventList
    *
    * @param runstart_i64 : absolute run start time in int64_t format with unit
    *nanosecond
    * @param pulsedt : pulse length in int64_t format with unit nanosecond
   */
  EventWorkspace_sptr createEventWorkspaceDirect(int64_t runstart_i64,
                                                 int64_t pulsedt) {
    // Create an EventWorkspace with 10 banks with 1 detector each.  No events
    // is generated
    EventWorkspace_sptr eventWS =
        WorkspaceCreationHelper::createEventWorkspaceWithFullInstrument(10, 1,
                                                                        true);

    // L1 = 10
    const auto &spectrumInfo = eventWS->spectrumInfo();
    double l1 = spectrumInfo.l1();

    Kernel::DateAndTime runstart(runstart_i64);

    EventList fakeevlist = fake_uniform_time_sns_data(runstart_i64, pulsedt);

    // Set properties: (1) run_start time; (2) Ei
    eventWS->mutableRun().addProperty("run_start", runstart.toISO8601String(),
                                      true);

    double shift = 2.E-4;
    double ei = (l1 * l1 * PhysicalConstants::NeutronMass) /
                (shift * shift * 2. * PhysicalConstants::meV);

    eventWS->mutableRun().addProperty<double>("Ei", ei, true);

    // Add neutrons
    for (size_t i = 0; i < eventWS->getNumberHistograms(); i++) {
      auto &elist = eventWS->getSpectrum(i);

      for (size_t ievent = 0; ievent < fakeevlist.getNumberEvents(); ++ievent) {
        TofEvent tofevent = fakeevlist.getEvent(ievent);
        elist.addEventQuickly(tofevent);
      } // FOR each pulse
    }   // For each bank

    // double constshift = l1 / sqrt(ei * 2. * PhysicalConstants::meV /
    //                           PhysicalConstants::NeutronMass);

    return eventWS;
  }

  //----------------------------------------------------------------------------------------------
  /** Create an EventWorkspace to mimic direct inelastic scattering insturment.
   * This workspace has
    * @param runstart_i64 : absolute run start time in int64_t format with unit
   * nanosecond
    * @param pulsedt : pulse length in int64_t format with unit nanosecond
    */
  EventWorkspace_sptr createEventWorkspaceInDirect(int64_t runstart_i64,
                                                   int64_t pulsedt) {
    // Create an EventWorkspace with 10 banks with 1 detector each.  No events
    // is generated
    EventWorkspace_sptr eventWS =
        WorkspaceCreationHelper::createEventWorkspaceWithFullInstrument(10, 1,
                                                                        true);
    // Add EFixed to each detector
    const ParameterMap &pmap = eventWS->constInstrumentParameters();

    for (size_t i = 0; i < 10; ++i) {

      const auto &spectrumInfo = eventWS->spectrumInfo();

      const auto &det = spectrumInfo.detector(i);
      Parameter_sptr par = pmap.getRecursive(&det, "Efixed");
      if (par) {
        // No need to set up E-Fix
        // double efix = par->value<double>();
        ;
      } else {

        eventWS->setEFixed(det.getID(), 2.08);
      }
    }

    // Add neutrons
    EventList fakeevlist = fake_uniform_time_sns_data(runstart_i64, pulsedt);
    for (size_t i = 0; i < eventWS->getNumberHistograms(); i++) {
      auto &elist = eventWS->getSpectrum(i);

      for (size_t ievent = 0; ievent < fakeevlist.getNumberEvents(); ++ievent) {
        TofEvent tofevent = fakeevlist.getEvent(ievent);
        elist.addEventQuickly(tofevent);
      } // FOR each pulse
    }   // For each bank

    return eventWS;
  }

  //----------------------------------------------------------------------------------------------
  /** Create an EventWorkspace as diffractometer
   * @brief createEventWorkspaceElastic
   * @param runstart_i64
   * @param pulsedt
   * @return
   */
  EventWorkspace_sptr createEventWorkspaceElastic(int64_t runstart_i64,
                                                  int64_t pulsedt) {
    // Create an EventWorkspace with 10 banks with 1 detector each.  No events
    // is generated
    EventWorkspace_sptr eventWS =
        WorkspaceCreationHelper::createEventWorkspaceWithFullInstrument(10, 1,
                                                                        true);

    Kernel::DateAndTime runstart(runstart_i64);

    // Create 1000 events
    EventList fakeevlist = fake_uniform_time_sns_data(runstart_i64, pulsedt);

    // Set properties: (1) run_start time; (2) Ei
    eventWS->mutableRun().addProperty("run_start", runstart.toISO8601String(),
                                      true);

    // Add neutrons
    for (size_t i = 0; i < eventWS->getNumberHistograms(); i++) {
      auto &elist = eventWS->getSpectrum(i);

      for (size_t ievent = 0; ievent < fakeevlist.getNumberEvents(); ++ievent) {
        TofEvent tofevent = fakeevlist.getEvent(ievent);
        elist.addEventQuickly(tofevent);
      } // FOR each pulse
    }   // For each bank

    return eventWS;
  }

  //----------------------------------------------------------------------------------------------
  /** Create a  Splitter for output
   *  Region:
   * 0: pulse 0: 0 ~ 3+
   * 1: pulse 0: 3+ ~ pulse 1: 9+
   * 2: from pulse 2: 0 ~ 6+
   * -1: from pulse 2: 6+ ~ 9+
    * @param runstart_i64 : absolute run start time in int64_t format with unit
   * nanosecond
    * @param pulsedt : pulse length in int64_t format with unit nanosecond
    * @param todft : time interval between 2 adjacent event in same pulse in
   * int64_t format of unit nanosecond
    * @param numpulses : number of pulses in the event workspace
   */
  SplittersWorkspace_sptr createSplitter(int64_t runstart_i64, int64_t pulsedt,
                                         int64_t tofdt) {
    SplittersWorkspace_sptr splitterws =
        boost::shared_ptr<SplittersWorkspace>(new SplittersWorkspace);

    // 1. Splitter 0: 0 ~ 3+ (first pulse)
    int64_t t0 = runstart_i64;
    int64_t t1 = t0 + tofdt * 3 + tofdt / 2;
    Kernel::SplittingInterval interval0(t0, t1, 0);
    splitterws->addSplitter(interval0);

    // 2. Splitter 1: 3+ ~ 9+ (second pulse)
    t0 = t1;
    t1 = runstart_i64 + pulsedt + tofdt * 9 + tofdt / 2;
    Kernel::SplittingInterval interval1(t0, t1, 1);
    splitterws->addSplitter(interval1);

    // 3. Splitter 2: from 3rd pulse, 0 ~ 6+
    for (size_t i = 2; i < 5; i++) {
      t0 = runstart_i64 + i * pulsedt;
      t1 = runstart_i64 + i * pulsedt + 6 * tofdt + tofdt / 2;
      Kernel::SplittingInterval interval2(t0, t1, 2);
      splitterws->addSplitter(interval2);
    }

    return splitterws;
  }

  //----------------------------------------------------------------------------------------------
  /** Create a  Splitter for output
   *  Region:
   * 0: pulse 0: 0 ~ 3+
   * 1: pulse 0: 3+ ~ pulse 1: 9+
   * 2: from pulse 2: 0 ~ 6+
   * -1: from pulse 2: 6+ ~ 9+
   * @brief createMatrixSplitter
   * @param runstart_i64 : absolute run start time in int64_t format with unit
   * nanosecond
   * @param pulsedt: pulse length in int64_t format with unit nanosecond
   * @param tofdt: time interval between 2 adjacent event in same pulse in
   * int64_t format of unit nanosecond
   * @return
   */
  API::MatrixWorkspace_sptr
  createMatrixSplitter(int64_t runstart_i64, int64_t pulsedt, int64_t tofdt) {
    // Create vectors for the splitters
    std::vector<int64_t> time_vec;
    std::vector<int> index_vec;

    time_vec.push_back(runstart_i64);

    // Splitter 0: 0 ~ 3+ (first pulse)
    int64_t t1 = runstart_i64 + tofdt * 3 + tofdt / 2;
    time_vec.push_back(t1);
    index_vec.push_back(0);

    // Splitter 1: 3+ ~ 9+ (second pulse)
    int64_t t2 = runstart_i64 + pulsedt + tofdt * 9 + tofdt / 2;
    time_vec.push_back(t2);
    index_vec.push_back(1);

    // Splitter 2 and so on: from 3rd pulse, 0 ~ 6+
    for (size_t i = 2; i < 5; i++) {
      int64_t newT = runstart_i64 + i * pulsedt + 6 * tofdt + tofdt / 2;
      time_vec.push_back(newT);
      index_vec.push_back(2);
    }

    // Create the workspace and set it
    size_t size_x = time_vec.size();
    size_t size_y = index_vec.size();
    TS_ASSERT(size_x - size_y == 1);

    MatrixWorkspace_sptr splitterws =
        boost::dynamic_pointer_cast<MatrixWorkspace>(
            WorkspaceFactory::Instance().create("Workspace2D", 1, size_x,
                                                size_y));

    for (size_t ix = 0; ix < size_x; ++ix)
      splitterws->mutableX(0)[ix] = static_cast<double>(time_vec[ix]);
    for (size_t iy = 0; iy < size_y; ++iy)
      splitterws->mutableY(0)[iy] = static_cast<double>(index_vec[iy]);

    return splitterws;
  }

  /** Create splitters in TableWorkspace for output which is exactly as the
   * Matrix splitters
   *  Region:
   * 0: pulse 0: 0 ~ 3+
   * 1: pulse 0: 3+ ~ pulse 1: 9+
   * 2: from pulse 2: 0 ~ 6+
   * -1: from pulse 2: 6+ ~ 9+
   * @brief createMatrixSplitter
   * @param runstart_i64 : absolute run start time in int64_t format with unit
   * nanosecond
   * @param pulsedt: pulse length in int64_t format with unit nanosecond
   * @param tofdt: time interval between 2 adjacent event in same pulse in
   * int64_t format of unit nanosecond
   * @return
   */
  DataObjects::TableWorkspace_sptr
  createTableSplitters(int64_t runstart_i64, int64_t pulsedt, int64_t tofdt) {
    // create table workspace
    DataObjects::TableWorkspace_sptr tablesplitter =
        boost::make_shared<DataObjects::TableWorkspace>();
    tablesplitter->addColumn("double", "start");
    tablesplitter->addColumn("double", "stop");
    tablesplitter->addColumn("str", "target");

    // generate row by row
    // Splitter 0: 0 ~ 3+ (first pulse)
    size_t row_index = 0;
    int64_t t1 = runstart_i64 + tofdt * 3 + tofdt / 2;
    std::string itarget = "A";
    tablesplitter->appendRow();
    tablesplitter->cell<double>(row_index, 0) =
        static_cast<double>(runstart_i64) * 1.0E-9;
    tablesplitter->cell<double>(row_index, 1) = static_cast<double>(t1) * 1.E-9;
    tablesplitter->cell<std::string>(row_index, 2) = itarget;

    // Splitter 1: 3+ ~ 9+ (second pulse)
    ++row_index;
    int64_t t2 = runstart_i64 + pulsedt + tofdt * 9 + tofdt / 2;
    itarget = "B";
    tablesplitter->appendRow();
    tablesplitter->cell<double>(row_index, 0) =
        static_cast<double>(t1) * 1.0E-9;
    tablesplitter->cell<double>(row_index, 1) = static_cast<double>(t2) * 1.E-9;
    tablesplitter->cell<std::string>(row_index, 2) = itarget;

    // Splitter 2 and so on: from 3rd pulse, 0 ~ 6+
    int64_t lastT = t2;
    for (size_t i = 2; i < 5; i++) {
      ++row_index;
      itarget = "C";
      int64_t newT = runstart_i64 + i * pulsedt + 6 * tofdt + tofdt / 2;
      tablesplitter->appendRow();
      tablesplitter->cell<double>(row_index, 0) =
          static_cast<double>(lastT) * 1.0E-9;
      tablesplitter->cell<double>(row_index, 1) =
          static_cast<double>(newT) * 1.E-9;
      tablesplitter->cell<std::string>(row_index, 2) = itarget;
      lastT = newT;
    }

    return tablesplitter;
  }

  //----------------------------------------------------------------------------------------------
  /** Create a Splitter for fast fequency log for output
    * The splitter is within every pulse.  2 groups of splitters are created.
    *In each pulse
    * 1. group 0: 0.2 dT ~ 0.4 dT    dT = pulsedt
    * 2. group 1: 0.6 dT ~ 0.8 dT
    *
    * @param runstart_i64 : absolute run start time in int64_t format with unit
    *nanosecond
    * @param pulsedt : pulse length in int64_t format with unit nanosecond
    * @param todft : time interval between 2 adjacent event in same pulse in
    *int64_t format of unit nanosecond
    * @param numpulses : number of pulses in the event workspace
   */
  SplittersWorkspace_sptr createFastFreqLogSplitter(int64_t runstart_i64,
                                                    int64_t pulsedt,
                                                    int64_t tofdt,
                                                    size_t numpulses) {

    UNUSED_ARG(tofdt);

    // 1. Create an empty splitter workspace
    SplittersWorkspace_sptr splitterws =
        boost::shared_ptr<SplittersWorkspace>(new SplittersWorkspace);

    // 2. Create splitters
    for (size_t i = 0; i < numpulses; ++i) {
      int64_t t0a = runstart_i64 + static_cast<int64_t>(i) * pulsedt +
                    static_cast<int64_t>(static_cast<double>(pulsedt) * 0.2);
      int64_t tfa = runstart_i64 + static_cast<int64_t>(i) * pulsedt +
                    static_cast<int64_t>(static_cast<double>(pulsedt) * 0.4);
      Kernel::SplittingInterval interval_a(t0a, tfa, 0);

      int64_t t0b = runstart_i64 + static_cast<int64_t>(i) * pulsedt +
                    static_cast<int64_t>(static_cast<double>(pulsedt) * 0.6);
      int64_t tfb = runstart_i64 + static_cast<int64_t>(i) * pulsedt +
                    static_cast<int64_t>(static_cast<double>(pulsedt) * 0.8);
      Kernel::SplittingInterval interval_b(t0b, tfb, 1);

      splitterws->addSplitter(interval_a);
      splitterws->addSplitter(interval_b);
    }

    return splitterws;
  }

  //----------------------------------------------------------------------------------------------
  /** Create the time correction table
    */
  TableWorkspace_sptr createTimeCorrectionTable(MatrixWorkspace_sptr inpws) {
    // 1. Generate an empty table
    auto corrtable = boost::make_shared<TableWorkspace>();
    corrtable->addColumn("int", "DetectorID");
    corrtable->addColumn("double", "Correction");

    // 2. Add rows
    const auto &detectorInfo = inpws->detectorInfo();
    const auto detids = detectorInfo.detectorIDs();
    for (size_t i = 0; i < detids.size(); ++i) {
      int detid = detids[i];
      double factor = 0.75;
      TableRow newrow = corrtable->appendRow();
      newrow << detid << factor;
    }

    return corrtable;
  }

  //----------------------------------------------------------------------------------------------
  /** Fake uniform time data more close to SNS case
    * A list of 1000 events
    * Pulse length: 1000000 * nano-second
    */
  EventList fake_uniform_time_sns_data(int64_t runstart, int64_t pulselength) {
    // Clear the list
    EventList el = EventList();

    // Create some mostly-reasonable fake data.
    unsigned seed1 = 1;
    std::minstd_rand0 g1(seed1);

    for (int time = 0; time < 1000; time++) {
      // All pulse times from 0 to 999 in seconds
      Kernel::DateAndTime pulsetime(
          static_cast<int64_t>(time * pulselength + runstart));
      double tof = static_cast<double>(g1() % 1000);
      el += TofEvent(tof, pulsetime);
    }

    return el;
  }

  /** Create a matrix splitters workspace for elastic correction
   * @brief createMatrixSplittersElastic
   * @return
   */
  API::MatrixWorkspace_sptr createMatrixSplittersElastic() {
    MatrixWorkspace_sptr spws = boost::dynamic_pointer_cast<MatrixWorkspace>(
        WorkspaceFactory::Instance().create("Workspace2D", 1, 11, 10));

    auto &vec_splitTimes = spws->mutableX(0);
    auto &vec_splitGroup = spws->mutableY(0);

    vec_splitTimes[0] = 1000000;
    vec_splitTimes[1] = 1300000;
    vec_splitTimes[2] = 2000000;
    vec_splitTimes[3] = 2190000;
    vec_splitTimes[4] = 4000000;
    vec_splitTimes[5] = 5000000;
    vec_splitTimes[6] = 5500000;
    vec_splitTimes[7] = 7000000;
    vec_splitTimes[8] = 8000000;
    vec_splitTimes[9] = 9000000;
    vec_splitTimes[10] = 10000000;

    vec_splitGroup[0] = 2;
    vec_splitGroup[1] = 5;
    vec_splitGroup[2] = 4;
    vec_splitGroup[3] = -1;
    vec_splitGroup[4] = 6;
    vec_splitGroup[5] = 7;
    vec_splitGroup[6] = 8;
    vec_splitGroup[7] = -1;
    vec_splitGroup[8] = 1;
    vec_splitGroup[9] = 3;

    return spws;
  }

  API::MatrixWorkspace_sptr createMatrixSplittersDG() {
    MatrixWorkspace_sptr spws = boost::dynamic_pointer_cast<MatrixWorkspace>(
        WorkspaceFactory::Instance().create("Workspace2D", 1, 11, 10));

    auto &vec_splitTimes = spws->mutableX(0);
    auto &vec_splitGroup = spws->mutableY(0);

    vec_splitTimes[0] = 1000000;
    vec_splitTimes[1] = 1300000; // Rule in  1,339,000
    vec_splitTimes[2] = 2000000;
    vec_splitTimes[3] = 2190000; // Rule out 2,155,000
    vec_splitTimes[4] = 4000000;
    vec_splitTimes[5] = 5000000;
    vec_splitTimes[6] = 5500000; // Rule in  5,741,000
    vec_splitTimes[7] = 7000000;
    vec_splitTimes[8] = 8000000;
    vec_splitTimes[9] = 9000000;
    vec_splitTimes[10] = 10000000;

    vec_splitGroup[0] = 2;
    vec_splitGroup[1] = 5;
    vec_splitGroup[2] = 4;
    vec_splitGroup[3] = -1;
    vec_splitGroup[4] = 6;
    vec_splitGroup[5] = 7;
    vec_splitGroup[6] = 8;
    vec_splitGroup[7] = -1;
    vec_splitGroup[8] = 1;
    vec_splitGroup[9] = 3;

    return spws;
  }
};

#endif /* MANTID_ALGORITHMS_FILTEREVENTSTEST_H_ */
