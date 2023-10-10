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
#pragma once

#include "gstlearn_export.hpp"

class ALinearOp;

class GSTLEARN_EXPORT LogStats {

public:
  LogStats(bool mustPrint = false);
  LogStats(const LogStats &m);
  LogStats& operator=(const LogStats &m);
  virtual ~LogStats();

  // The four following methods are considered as const not to destroy constness of calling functions
  void incrementStatsDirect(double time) const;
  void incrementStatsInverseCG(int niter, double time) const;
  void incrementStatsInverseChol(double time) const;
  void incrementStatsInversePoly(double time) const;
  void incrementStatsSimulate(double time) const;
  void incrementStatsCholesky(double time) const;

  void statsShow(void) const;

  void mustShowStats(bool mustPrint) const { _mustPrint = mustPrint; }
  bool isMustPrint() const { return _mustPrint; }

private:
  mutable bool       _mustPrint;

  mutable int        _directNumber;
  mutable double     _directTime;

  mutable int        _inverseCGNIter;
  mutable int        _inverseCGNumber;
  mutable double     _inverseCGTime;

  mutable int        _inverseCholNumber;
  mutable double     _inverseCholTime;

  mutable int        _inversePolyNumber;
  mutable double     _inversePolyTime;

  mutable int        _simulateNumber;
  mutable double     _simulateTime;

  mutable int        _choleskyNumber;
  mutable double     _choleskyTime;
};
