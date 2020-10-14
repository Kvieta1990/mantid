// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2012 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

/*
    This file contains the DLLExport/DLLImport linkage configuration for the
    MDAlgorithms library

    @author Martyn Gigg, Tessella plc
*/
#include "MantidKernel/System.h"

#ifdef IN_MANTID_MDALGORITHMS
#define MANTID_MDALGORITHMS_DLL DLLExport
#define EXTERN_MANTID_MDALGORITHMS
#else
#define MANTID_MDALGORITHMS_DLL DLLImport
#define EXTERN_MANTID_MDALGORITHMS EXTERN_IMPORT
#endif /* IN_MANTID_MDALGORITHMS*/
