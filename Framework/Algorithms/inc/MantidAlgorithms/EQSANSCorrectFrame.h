// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAPI/Algorithm.h"
#include "MantidAPI/MatrixWorkspace_fwd.h"
#include "MantidAlgorithms/DllConfig.h"
#include "MantidDataObjects/EventWorkspace.h"

namespace Mantid {
namespace Algorithms {
/**
    Apply correction to EQSANS data to account for its TOF structure. The
   algorithm modifies the TOF values to correct for the fact that T_0 is not
   properly recorded by the DAS.

    File change history is stored at: <https://github.com/mantidproject/mantid>
    Code Documentation is available at: <http://doxygen.mantidproject.org>
*/

class MANTID_ALGORITHMS_DLL EQSANSCorrectFrame : public API::Algorithm {
public:
  /// Default constructor
  EQSANSCorrectFrame();
  /// Algorithm's name
  const std::string name() const override { return "EQSANSCorrectFrame"; }
  /// Summary of algorithms purpose
  const std::string summary() const override {
    return "Corrects the TOF of raw EQSANS data. This algorithm needs to be "
           "run once on every data set.";
  }

  /// Algorithm's version
  int version() const override { return (1); }
  /// Algorithm's category for identification
  const std::string category() const override { return "SANS"; }

private:
  /// Initialisation code
  void init() override;
  /// Execution code
  void exec() override;
};

} // namespace Algorithms
} // namespace Mantid
