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
#include "Basic/Timer.hpp"
#include "Basic/AStringable.hpp"

Timer::Timer()
{
  reset();
}

Timer::Timer(const Timer &m)
    : _refTime(m._refTime)
{

}

Timer& Timer::operator=(const Timer &m)
{
  if (this != &m)
  {
    _refTime = m._refTime;
  }
  return *this;
}

Timer::~Timer()
{

}

/**
 * Defines the reference Timer (beginning of a Timer chunk)
 */
void Timer::reset()
{
  _refTime = hrc::now();
}

/**
 * Displays the time elapsed (in s) since the reference Timer
 * @param title Title used for the internal display
 * @param flag_reset True if the Reference must be set to current Time
 */
void Timer::displayIntervalSeconds(const String& title, bool flag_reset)
{
  double sec = getIntervalSeconds(flag_reset);
  displaySeconds(title, sec);
}

double Timer::getIntervalSeconds(bool flag_reset)
{
  auto newTime = hrc::now();
  sec fs = newTime - _refTime;
  if (flag_reset) _refTime = newTime;
  return fs.count();
}

void Timer::displaySeconds(const String& title, double sec)
{
  if (! title.empty())
    message("%s: %6.2lf s.\n",title.c_str(),sec);
  else
    message("Timer: %6.2lf s.\n",sec);

}

/**
 * Displays the time elapsed (in ms) since the reference Timer
 * @param title Title used for the internal display
 * @param flag_reset True if the Reference must be set to current Time
 */
void Timer::displayIntervalMilliseconds(const String& title, bool flag_reset)
{
  double msec = getIntervalMilliseconds(flag_reset);
  displayMilliseconds(title, msec);
}

double Timer::getIntervalMilliseconds(bool flag_reset)
{
  auto newTime = hrc::now();
  sec fs = newTime - _refTime;
  ms inter = std::chrono::duration_cast<ms>(fs);
  if (flag_reset) _refTime = newTime;
  return ((double) inter.count());
}

void Timer::displayMilliseconds(const String& title, double msec)
{
  if (! title.empty())
    message("%s: %6.2lf ms.\n",title.c_str(),msec);
  else
    message("Timer: %6.2lf ms.\n",msec);
}
