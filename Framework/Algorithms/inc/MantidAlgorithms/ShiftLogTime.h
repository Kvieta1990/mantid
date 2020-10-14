// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2011 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/Algorithm.h"
#include "MantidAlgorithms/DllConfig.h"

namespace Mantid {
namespace Algorithms {

/** ShiftLogTime : TODO: DESCRIPTION

  @date 2011-08-26
*/
class MANTID_ALGORITHMS_DLL ShiftLogTime : public API::Algorithm {
public:
  const std::string name() const override;
  int version() const override;
  const std::vector<std::string> seeAlso() const override {
    return {"CreateLogTimeCorrection", "ChangePulsetime", "ChangeLogTime"};
  }
  const std::string category() const override;

  /// Summary of algorithms purpose
  const std::string summary() const override {
    return "Shifts the indexes of the specified log. This will make the log "
           "shorter by the specified shift.";
  }

private:
  /// Sets documentation strings for this algorithm

  /// Initialise the properties
  void init() override;
  /// Run the algorithm
  void exec() override;
};

} // namespace Algorithms
} // namespace Mantid
