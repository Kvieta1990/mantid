// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidKernel/DllConfig.h"
#include "MantidKernel/IPropertySettings.h"

namespace Mantid {
namespace Kernel {

/** This property setting object makes a property invisible in the GUI.
 */
class MANTID_KERNEL_DLL InvisibleProperty : public IPropertySettings {
public:
  bool isVisible(const IPropertyManager *) const override;
  IPropertySettings *clone() const override;
};

} // namespace Kernel
} // namespace Mantid
