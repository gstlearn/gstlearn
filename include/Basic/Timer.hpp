/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
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

  void displayIntervalSeconds(const String& title = String(), bool flag_reset = true);
  double getIntervalSeconds(bool flag_reset = true);
  void displaySeconds(const String& title, double sec);

  void displayIntervalMilliseconds(const String& title = String(), bool flag_reset = true);
  double getIntervalMilliseconds(bool flag_reset = true);
  void displayMilliseconds(const String& title, double msec);

private:
  hrc::time_point _refTime;
};
