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
#include "Basic/AStringable.hpp"

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
 * @param flag_reset True if the Reference must be set to current Time
 * @param verbose True if a message should be issued
 * @param title Title used for the internal display
 * @return Timer elapsed (ms)
 */
double Timer::getTimerInterval(bool flag_reset,
                               bool verbose,
                               const String& title)
{
  clock_t newTimer  = clock();
  clock_t diffTimer = newTimer - _refTimer;
  if (flag_reset) _refTimer = newTimer;
  double msec = (double) (diffTimer) / CLOCKS_PER_SEC;

  if (verbose)
  {
    if (! title.empty())
      message("%s: %6.2lf ms.\n",title.c_str(),msec);
    else
      message("Timer: %6.2lf ms.\n",msec);
  }

  return msec;
}
