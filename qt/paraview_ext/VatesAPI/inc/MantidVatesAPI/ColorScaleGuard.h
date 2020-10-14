// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once
#include "MantidKernel/Logger.h"
#include "MantidKernel/System.h"

namespace Mantid {
namespace VATES {

class DLLExport ColorScaleLock {
public:
  bool isLocked() { return m_isLocked; }
  void lock() { m_isLocked = true; }
  void unlock() { m_isLocked = false; }

private:
  bool m_isLocked{false};
};

class DLLExport ColorScaleLockGuard {
public:
  ColorScaleLockGuard(ColorScaleLock *lock) {
    if (!lock || lock->isLocked()) {
      m_lock = nullptr;
    } else {
      m_lock = lock;
      m_lock->lock();
    }
  }

  ~ColorScaleLockGuard() {
    if (m_lock) {
      m_lock->unlock();
    }
  }

private:
  ColorScaleLock *m_lock;
};
} // namespace VATES
} // namespace Mantid