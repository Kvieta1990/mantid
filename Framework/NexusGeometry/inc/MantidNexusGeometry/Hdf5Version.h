// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidNexusGeometry/DllConfig.h"

#include <cstdint>

namespace Mantid {
namespace NexusGeometry {
namespace Hdf5Version {

// Create a comparable version number as a single integer
uint32_t makeHdf5VersionNumber(uint32_t maj, uint32_t min, uint32_t relnum);

// Check if current version of hdf5 supports variable length strings
MANTID_NEXUSGEOMETRY_DLL bool checkVariableLengthStringSupport();

} // namespace Hdf5Version
} // namespace NexusGeometry
} // namespace Mantid
