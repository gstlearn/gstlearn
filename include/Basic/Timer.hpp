/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include <time.h>

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
  void Interval(const String& title = String(), bool flag_reset = true);
  double getInterval(bool flag_reset = true);
  void display(const String& title, double msec);

private:
  clock_t _refTimer;
};
