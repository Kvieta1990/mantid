// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/IFunctionMD.h"
#include <cxxtest/TestSuite.h>

using namespace Mantid;
using namespace Mantid::API;

class IFunctionMDTest : public CxxTest::TestSuite {
public:
  void testIFunction() {
    // What does this test do??
    int i = 1;
    TS_ASSERT(i);

    // const int nx = 10;
    // const int ny = 10;
    // Mantid::DataObjects::Workspace2D_sptr ws = Create2DWorkspace(nx,ny);

    // for(int i=0;i<nx;++i)
    // {
    //   for(int j=0;j<ny;++j)
    //   {
    //   }
    // }

    // std::cerr<<"\nn="<<ws->axes()<<"\n";
  }
};
