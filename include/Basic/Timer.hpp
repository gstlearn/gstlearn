/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include <chrono>

typedef std::chrono::high_resolution_clock hrc;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<double> sec;

/**
 * Allow calculation of Timer spent in a portion of the code
 */
class GSTLEARN_EXPORT Timer
{
public:
  Timer();
  Timer(const Timer &m);
  Timer& operator= (const Timer &m);
  virtual ~Timer();

  void reset();

  void displayIntervalSeconds(const String& title = String(),
                              int expected_time = -1,
                              bool flag_reset = true);
  double getIntervalSeconds(bool flag_reset = true);
  void displaySeconds(const String& title,
                      double sec,
                      int expected_time = -1);

  void displayIntervalMilliseconds(const String& title = String(),
                                   int expected_time = -1,
                                   bool flag_reset = true);
  double getIntervalMilliseconds(bool flag_reset = true);
  void displayMilliseconds(const String& title,
                           double msec,
                           int expected_time = -1);

private:
  hrc::time_point _refTime;
};
