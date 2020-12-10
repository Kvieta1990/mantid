// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/Run.h"
#include "MantidKernel/System.h"
#include "MantidKernel/Timer.h"
#include <cxxtest/TestSuite.h>

#include "MantidDataHandling/LoadDetectorsGroupingFile.h"
#include "MantidDataHandling/SaveDetectorsGrouping.h"
#include "Poco/File.h"

using namespace Mantid;
using namespace Mantid::API;
using namespace Mantid::DataHandling;
using namespace Mantid::DataObjects;

namespace {

const std::string TEMP_SAVED_MAP_WS = "test_grouping_workspace";
const std::string TEMP_LOADED_MAP_WS = "test_grouping_loaded";
const std::string TEMP_MAP_FILE_NAME = "temp_grouping.map";

void createGroupingWorkspace(const std::string &instrumentName,
                             const std::string &analyser,
                             const std::string &customGrouping) {
  auto creator = AlgorithmManager::Instance().create("CreateGroupingWorkspace");
  creator->initialize();
  creator->setProperty("InstrumentName", instrumentName);
  creator->setProperty("ComponentName", analyser);
  creator->setProperty("CustomGroupingString", customGrouping);
  creator->setProperty("OutputWorkspace", TEMP_SAVED_MAP_WS);
  creator->execute();
}

ITableWorkspace_sptr compareWorkspaces(const MatrixWorkspace_sptr &workspace1,
                                       const MatrixWorkspace_sptr &workspace2,
                                       double tolerance = 0.000001) {
  auto compareAlg = AlgorithmManager::Instance().create("CompareWorkspaces");
  compareAlg->setProperty("Workspace1", workspace1);
  compareAlg->setProperty("Workspace2", workspace2);
  compareAlg->setProperty("Tolerance", tolerance);
  compareAlg->execute();
  return compareAlg->getProperty("Messages");
}

} // namespace

class SaveDetectorsGroupingTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static SaveDetectorsGroupingTest *createSuite() {
    return new SaveDetectorsGroupingTest();
  }
  static void destroySuite(SaveDetectorsGroupingTest *suite) { delete suite; }

  void test_Initialize() {
    SaveDetectorsGrouping savegroup;
    savegroup.initialize();
    TS_ASSERT(savegroup.isInitialized());
  }

  void test_SaveXMLFile() {

    // 1. Get an object for SaveDetectorsGrouping
    SaveDetectorsGrouping savegroup;
    savegroup.initialize();

    // 2. Create a Workspace from an XML file
    LoadDetectorsGroupingFile load;
    load.initialize();

    TS_ASSERT(load.setProperty("InputFile", "vulcangroup.xml"));
    TS_ASSERT(load.setProperty("OutputWorkspace", "Vulcan_Group"));

    load.execute();
    TS_ASSERT(load.isExecuted());

    DataObjects::GroupingWorkspace_sptr gws =
        std::dynamic_pointer_cast<DataObjects::GroupingWorkspace>(
            API::AnalysisDataService::Instance().retrieve("Vulcan_Group"));

    // 3. Save
    savegroup.setProperty("InputWorkspace", gws);
    savegroup.setProperty("OutputFile", "grouptemp.xml");

    savegroup.execute();
    TS_ASSERT(savegroup.isExecuted());

    // Retrieve the full path to file
    std::string file1 = savegroup.getPropertyValue("OutputFile");

    // 4. Verify
    LoadDetectorsGroupingFile load2;
    load2.initialize();

    TS_ASSERT(load2.setProperty("InputFile", file1));
    TS_ASSERT(load2.setProperty("OutputWorkspace", "Vulcan_Group2"));

    load2.execute();
    TS_ASSERT(load2.isExecuted());

    DataObjects::GroupingWorkspace_sptr gws2 =
        std::dynamic_pointer_cast<DataObjects::GroupingWorkspace>(
            API::AnalysisDataService::Instance().retrieve("Vulcan_Group2"));

    TS_ASSERT_DELTA(gws2->y(0)[0], 1.0, 1.0E-5);
    // TS_ASSERT_DELTA(gws2->y(3695)[0], 2.0, 1.0E-5);
    TS_ASSERT_DELTA(gws2->y(3696)[0], 2.0, 1.0E-5);
    // TS_ASSERT_DELTA(gws2->y(7000)[0], 0.0, 1.0E-5);

    // 5. Clear
    Poco::File file(file1);
    file.remove(false);

    API::AnalysisDataService::Instance().remove("Vulcan_Group");
    API::AnalysisDataService::Instance().remove("Vulcan_Group2");
  }

  void test_saving_a_MAP_file() {
    auto &ads = AnalysisDataService::Instance();

    createGroupingWorkspace("IRIS", "graphite", "3-5,6+7,8:10");

    TS_ASSERT(ads.doesExist(TEMP_SAVED_MAP_WS));
    const auto groupingToSave =
        ads.retrieveWS<GroupingWorkspace>(TEMP_SAVED_MAP_WS);

    SaveDetectorsGrouping saver;
    saver.initialize();
    saver.setProperty("InputWorkspace", groupingToSave);
    saver.setProperty("OutputFile", TEMP_MAP_FILE_NAME);
    saver.execute();
    TS_ASSERT(saver.isExecuted());

    LoadDetectorsGroupingFile loader;
    loader.initialize();
    loader.setProperty("InputFile", TEMP_MAP_FILE_NAME);
    loader.setProperty("OutputWorkspace", TEMP_LOADED_MAP_WS);
    loader.execute();
    TS_ASSERT(loader.isExecuted());

    const auto groupingLoaded =
        ads.retrieveWS<GroupingWorkspace>(TEMP_LOADED_MAP_WS);

    TS_ASSERT(!compareWorkspaces(groupingToSave, groupingLoaded));

    Poco::File file(TEMP_MAP_FILE_NAME);
    file.remove(false);
    ads.remove(TEMP_SAVED_MAP_WS);
    ads.remove(TEMP_LOADED_MAP_WS);
  }

  void test_SaveNamingAndDescription() {
    std::string testWs = "GroupingWorkspace";

    // Load the grouping to test with
    LoadDetectorsGroupingFile load;
    load.initialize();

    TS_ASSERT(load.setProperty("InputFile", "MUSRGrouping.xml"));
    TS_ASSERT(load.setProperty("OutputWorkspace", testWs));

    load.execute();
    TS_ASSERT(load.isExecuted());

    // Save the test workspace
    SaveDetectorsGrouping save;
    save.initialize();

    TS_ASSERT(save.setProperty("InputWorkspace", testWs));
    TS_ASSERT(save.setProperty("OutputFile", "grouptemp.xml"));

    save.execute();
    TS_ASSERT(save.isExecuted());

    // Get full path of the file
    std::string testFile = save.getPropertyValue("OutputFile");

    // Remove saved workspace
    API::AnalysisDataService::Instance().remove(testWs);

    // Load it again to verify
    LoadDetectorsGroupingFile load2;
    load2.initialize();

    TS_ASSERT(load2.setProperty("InputFile", testFile));
    TS_ASSERT(load2.setProperty("OutputWorkspace", testWs));

    load2.execute();
    TS_ASSERT(load2.isExecuted());

    // Get GroupingWorkspace instance
    auto gws = std::dynamic_pointer_cast<DataObjects::GroupingWorkspace>(
        API::AnalysisDataService::Instance().retrieve(testWs));

    // Check that description was saved
    TS_ASSERT_EQUALS(gws->run().getProperty("Description")->value(),
                     "musr longitudinal (64 detectors)");

    // Check that group names were saved
    TS_ASSERT_EQUALS(gws->run().getProperty("GroupName_1")->value(), "fwd");
    TS_ASSERT_EQUALS(gws->run().getProperty("GroupName_2")->value(), "bwd");

    // Clean-up
    API::AnalysisDataService::Instance().remove(testWs);

    Poco::File file(testFile);
    file.remove();
  }
};
