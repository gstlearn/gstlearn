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

#include <time.h>

/**
 * Allow calculation of Timer spent in a portion of the code
 */
class Timer
{
public:
  Timer();
  Timer(const Timer &m);
  Timer& operator= (const Timer &m);
  virtual ~Timer();

  void reset();
  int  getTimerInterval(bool flag_reset = false);

private:
  clock_t _refTimer;
};
