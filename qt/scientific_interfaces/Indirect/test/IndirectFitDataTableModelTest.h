// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2019 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "IndirectFitAnalysisTab.h"
#include "MantidAPI/FunctionFactory.h"
#include <cxxtest/TestSuite.h>
#include <gmock/gmock.h>

#include "IndirectFitDataTableModel.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidTestHelpers/IndirectFitDataCreationHelper.h"

using namespace MantidQt::CustomInterfaces::IDA;
using namespace MantidQt::MantidWidgets;

class IndirectFitDataTableModelTest : public CxxTest::TestSuite {
public:
  IndirectFitDataTableModelTest() = default;

  void setUp() override {
    m_fitData = std::make_unique<IndirectFitDataTableModel>();
    auto resolutionWorkspace =
        Mantid::IndirectFitDataCreationHelper::createWorkspace(4, 5);
    auto dataWorkspace =
        Mantid::IndirectFitDataCreationHelper::createWorkspace(4, 5);
    Mantid::API::AnalysisDataService::Instance().addOrReplace(
        "resolution workspace", std::move(resolutionWorkspace));
    Mantid::API::AnalysisDataService::Instance().addOrReplace(
        "data workspace", std::move(dataWorkspace));
    m_fitData->addWorkspace("data workspace");
    m_fitData->setResolution("resolution workspace", TableDatasetIndex{0});
  }

  void tearDown() override { AnalysisDataService::Instance().clear(); }

  void test_hasWorkspace_returns_true_for_ws_in_model() {
    TS_ASSERT(m_fitData->hasWorkspace("data workspace"));
  }

  void test_hasWorkspace_returns_false_for_ws_in_model() {
    TS_ASSERT(!m_fitData->hasWorkspace("fake workspace"));
  }

  void test_getWorkspace_returns_nullptr_is_outside_of_range() {
    TS_ASSERT_EQUALS(m_fitData->getWorkspace(TableDatasetIndex{1}), nullptr);
  }

  void test_getWorkspace_returns_ws_in_range() {
    TS_ASSERT_EQUALS(m_fitData->getWorkspace(TableDatasetIndex{0})->getName(),
                     "data workspace");
  }

  void test_getSpectra_returns_empty_spectra_is_outside_of_range() {
    TS_ASSERT_EQUALS(m_fitData->getSpectra(TableDatasetIndex{1}).getString(),
                     "");
  }

  void test_getSpectra_returns_spectra_in_range() {
    TS_ASSERT_EQUALS(m_fitData->getSpectra(TableDatasetIndex{0}).getString(),
                     "0-3");
  }

  void test_isMultiFit_returns_false_for_more_single_ws() {
    TS_ASSERT(!m_fitData->isMultiFit());
  }

  void test_isMultiFit_returns_true_for_more_than_one_ws() {
    auto dataWorkspace =
        Mantid::IndirectFitDataCreationHelper::createWorkspace(4, 5);
    Mantid::API::AnalysisDataService::Instance().addOrReplace(
        "data workspace 2", std::move(dataWorkspace));
    m_fitData->addWorkspace("data workspace 2");
    TS_ASSERT(m_fitData->isMultiFit());
  }

  void test_getNumberOfWorkspaces_returns_correct_number_of_workspaces() {
    TS_ASSERT_EQUALS(m_fitData->getNumberOfWorkspaces(), 1);
    auto dataWorkspace =
        Mantid::IndirectFitDataCreationHelper::createWorkspace(4, 5);
    Mantid::API::AnalysisDataService::Instance().addOrReplace(
        "data workspace 2", std::move(dataWorkspace));
    m_fitData->addWorkspace("data workspace 2");
    TS_ASSERT_EQUALS(m_fitData->getNumberOfWorkspaces(), 2);

  }

  void test_getNumberOfSpectra_returns_correct_number_of_spectra() {
    TS_ASSERT_EQUALS(m_fitData->getNumberOfSpectra(TableDatasetIndex{0}), 4);
    auto dataWorkspace =
        Mantid::IndirectFitDataCreationHelper::createWorkspace(5, 5);
    Mantid::API::AnalysisDataService::Instance().addOrReplace(
        "data workspace 2", std::move(dataWorkspace));
    m_fitData->addWorkspace("data workspace 2");
    TS_ASSERT_EQUALS(m_fitData->getNumberOfSpectra(TableDatasetIndex{1}), 5);
  }

  void test_getNumberOfSpectra_raise_error_when_out_of_ws_range() {
    TS_ASSERT_EQUALS(m_fitData->getNumberOfSpectra(TableDatasetIndex{0}), 4);
    TS_ASSERT_THROWS(m_fitData->getNumberOfSpectra(TableDatasetIndex{1}),
                     const std::runtime_error &);
  }

  void test_that_getResolutionsForFit_return_correctly() {
    auto resolutionVector = m_fitData->getResolutionsForFit();

    TS_ASSERT_EQUALS(resolutionVector[2].first, "resolution workspace");
    TS_ASSERT_EQUALS(resolutionVector[2].second, 2);
  }

  void
  test_that_getResolutionsForFit_return_correctly_if_resolution_workspace_removed() {
    Mantid::API::AnalysisDataService::Instance().clear();

    auto resolutionVector = m_fitData->getResolutionsForFit();

    TS_ASSERT_EQUALS(resolutionVector[2].first, "");
    TS_ASSERT_EQUALS(resolutionVector[2].second, 2);
  }

  void test_can_set_spectra_on_existing_workspace() {
    m_fitData->setSpectra("1", TableDatasetIndex{0});

    TS_ASSERT_EQUALS(m_fitData->getSpectra(TableDatasetIndex{0}),
                     FunctionModelSpectra("1"));
  }

  void test_that_setting_spectra_on_non_existent_workspace_throws_exception() {
    TS_ASSERT_THROWS(m_fitData->setSpectra("1", TableDatasetIndex{1}),
                     const std::out_of_range &)
    TS_ASSERT_THROWS(
        m_fitData->setSpectra(FunctionModelSpectra("1"), TableDatasetIndex{1}),
        const std::out_of_range &)
  }

  void test_that_setting_startX_on_non_existent_workspace_throws_exception() {
    TS_ASSERT_THROWS(m_fitData->setStartX(0, TableDatasetIndex{1}),
                     const std::out_of_range &)
    TS_ASSERT_THROWS(
        m_fitData->setStartX(0, TableDatasetIndex{1},
                             MantidQt::MantidWidgets::WorkspaceIndex{10}),
        const std::out_of_range &)
  }

private:
  std::unique_ptr<IIndirectFitDataTableModel> m_fitData;
};
