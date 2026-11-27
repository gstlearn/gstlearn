/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2026) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <string>
#include <unordered_map>

#ifndef GSTLEARN_DB_OPTIM
#define GSTLEARN_DB_OPTIM 0
#endif

/**
 * @brief This class provides a mechanism to count the number of times specific loops
 * or nodes are executed in the code.
 * It is designed to help identify performance bottlenecks and optimize critical sections of the code.
 * It is conditioned on the GSTLEARN_DB_OPTIM macro, which can be set to 0 (disabled) or 1 (enabled).
 */
namespace gstlrn
{
  inline constexpr bool areCountersActive = (GSTLEARN_DB_OPTIM != 0);

  /**
   * Loop / node pass counters (active when GSTLEARN_DB_OPTIM != 0)
   */
  class GSTLEARN_EXPORT DbOptim
  {
  public:
    static void incrementCounter(const String& name);
    static void resetCounters();
    static void displayCounters();

  private:
    static std::unordered_map<std::string, long long> _counters;
    static bool _atExitRegistered;
  };
} // namespace gstlrn
