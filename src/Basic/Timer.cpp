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
#include "Basic/Timer.hpp"

Timer::Timer()
{
  reset();
}

Timer::Timer(const Timer &m)
    : _refTimer(m._refTimer)
{

}

Timer& Timer::operator=(const Timer &m)
{
  if (this != &m)
  {
    _refTimer = m._refTimer;
  }
  return *this;
}

Timer::~Timer()
{

}

/**
 * Defines the reference Timer (begining of a Timer chunk)
 */
void Timer::reset()
{
  _refTimer = clock();
}

/**
 * Returns Timer elapsed (in ms) since the reference Timer
 * @return Timer elapsed (ms)
 */
int Timer::getTimerInterval(bool flag_reset)
{
  clock_t newTimer  = clock();
  clock_t diffTimer = newTimer - _refTimer;
  if (flag_reset) _refTimer = newTimer;
  int msec = diffTimer / CLOCKS_PER_SEC;
  return msec;
}
