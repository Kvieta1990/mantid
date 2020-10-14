// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidKernel/Timer.h"
#include <cxxtest/TestSuite.h>

#include "MantidAlgorithms/Pause.h"

using namespace Mantid;
using namespace Mantid::Algorithms;
using namespace Mantid::API;

class PauseTest : public CxxTest::TestSuite {
public:
  /** Not much to test, just that it runs */
  void test_exec() {
    Pause alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize())
    TS_ASSERT(alg.isInitialized())
    TS_ASSERT_THROWS_NOTHING(alg.setProperty("Duration", 0.1));
    TS_ASSERT_THROWS_NOTHING(alg.execute(););
    TS_ASSERT(alg.isExecuted());
  }
};
