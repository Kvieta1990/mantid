// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2012 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidKernel/System.h"

class vtkStructuredGrid;

namespace Mantid {
namespace VATES {

/** Generates a vtkStructuredGrid with a single point. Note that this is not a
 Null Object for a vtkDataSet.

 @date 15/07/2015
*/

class DLLExport vtkNullStructuredGrid {

public:
  vtkNullStructuredGrid();

  ~vtkNullStructuredGrid();

  vtkStructuredGrid *createNullData();
};
} // namespace VATES
} // namespace Mantid