// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidKernel/Matrix.h"

// Auxiliary class to store quadratic coefficients
// See J. Skilling and R. K. Bryan: "Maximum entropy image reconstruction:
// general algorithm" (1984)

struct QuadraticCoefficients {
  Mantid::Kernel::DblMatrix s1; // Quadratic coefficient S_mu
  Mantid::Kernel::DblMatrix c1; // Quadratic coefficient C_mu
  Mantid::Kernel::DblMatrix s2; // Quadratic coefficient g_mu_nu
  Mantid::Kernel::DblMatrix c2; // Quadratic coefficient M_mu_nu
};
