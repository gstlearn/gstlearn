/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "LinearOp/LogStats.hpp"
#include "Basic/AStringable.hpp"

LogStats::LogStats(bool mustPrint)
    : _mustPrint(mustPrint),
      _directNumber(0),
      _directTime(0.),
      _inverseCGNIter(0),
      _inverseCGNumber(0),
      _inverseCGTime(0.),
      _inverseCholNumber(0),
      _inverseCholTime(0.),
      _inversePolyNumber(0),
      _inversePolyTime(0.),
      _simulateNumber(0),
      _simulateTime(0.),
      _choleskyNumber(0),
      _choleskyTime(0.)
{
}

LogStats::LogStats(const LogStats &m)
    : _mustPrint(m._mustPrint),
      _directNumber(m._directNumber),
      _directTime(m._directTime),
      _inverseCGNIter(m._inverseCGNIter),
      _inverseCGNumber(m._inverseCGNumber),
      _inverseCGTime(m._inverseCGTime),
      _inverseCholNumber(m._inverseCholNumber),
      _inverseCholTime(m._inverseCholTime),
      _inversePolyNumber(m._inversePolyNumber),
      _inversePolyTime(m._inversePolyTime),
      _simulateNumber(m._simulateNumber),
      _simulateTime(m._simulateTime),
      _choleskyNumber(m._choleskyNumber),
      _choleskyTime(m._choleskyTime)
{
}

LogStats& LogStats::operator=(const LogStats &m)
{
  if (this != &m)
  {
    _mustPrint = m ._mustPrint;
    _directNumber = m._directNumber;
    _directTime = m._directTime;
    _inverseCGNIter = m._inverseCGNIter;
    _inverseCGNumber = m._inverseCGNumber;
    _inverseCGTime = m._inverseCGTime;
    _inverseCholNumber = m._inverseCholNumber;
    _inverseCholTime = m._inverseCholTime;
    _inversePolyNumber = m._inversePolyNumber;
    _inversePolyTime = m._inversePolyTime;
    _simulateNumber = m._simulateNumber;
    _simulateTime = m._simulateTime;
    _choleskyNumber = m._choleskyNumber;
    _choleskyTime = m._choleskyTime;
  }
  return *this;
}

LogStats::~LogStats()
{
}

void LogStats::incrementStatsDirect(double time) const
{
  _directNumber ++;
  _directTime   += time * 1000;
}

void LogStats::incrementStatsInverseCG(int niter, double time) const
{
  _inverseCGNumber ++;
  _inverseCGNIter  += niter;
  _inverseCGTime   += time * 1000;
}

void LogStats::incrementStatsInverseChol(double time) const
{
  _inverseCholNumber ++;
  _inverseCholTime   += time * 1000;
}

void LogStats::incrementStatsInversePoly(double time) const
{
  _inversePolyNumber ++;
  _inversePolyTime   += time * 1000;
}

void LogStats::incrementStatsSimulate(double time) const
{
  _simulateNumber ++;
  _simulateTime   += time * 1000;
}

void LogStats::incrementStatsCholesky(double time) const
{
  _choleskyNumber ++;
  _choleskyTime   += time * 1000;
}

/*****************************************************************************/
/*!
 **  Trigger the printout of the statistics
 **
 *****************************************************************************/
void LogStats::statsShow(void) const
{
  if (! _mustPrint) return;
  if (_choleskyNumber > 0)
    message("Statistics - Cholesky decomposition (*%d) : %lf ms\n",
            _choleskyNumber, _choleskyTime);
  if (_directNumber > 0)
    message("Statistics - Direct evaluation (*%d) : %lf ms\n",
            _directNumber, _directTime);
  if (_inverseCGNumber > 0)
    message("Statistics - Inverse using Conjugate Gradient (*%d) : %lf ms (%d iterations)\n",
            _inverseCGNumber, _inverseCGTime, _inverseCGNIter);
  if (_inverseCholNumber > 0)
    message("Statistics - Inverse using Cholesky decomposition (*%d) : %lf ms\n",
            _inverseCholNumber, _inverseCholTime);
  if (_inversePolyNumber > 0)
    message("Statistics - Inverse using Polynomial expansion (*%d) : %lf ms\n",
            _inversePolyNumber, _inversePolyTime);
  if (_simulateNumber > 0)
    message("Statistics - Simulation evaluation (*%d) : %lf ms\n",
            _simulateNumber, _simulateTime);
}

