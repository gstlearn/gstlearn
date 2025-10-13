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

#include "LinearOp/LogStats.hpp"
#include <vector>
#include "geoslib_define.h"

namespace gstlrn
{
  class LogStats;

class GSTLEARN_EXPORT ALinearOpMulti {

public:
  ALinearOpMulti(Id nitermax = 1000, double eps = EPSILON8);
  ALinearOpMulti(const ALinearOpMulti &m);
  ALinearOpMulti& operator=(const ALinearOpMulti &m);
  virtual ~ALinearOpMulti();

  void initLk(const std::vector<std::vector<double>> &inv, std::vector<std::vector<double>> &outv) const;
  virtual Id sizes() const = 0;
  virtual Id size(Id) const = 0;

  void setNIterMax(Id nitermax) { _nIterMax = nitermax; }
  void setNIterRestart(Id niterrestart) { _nIterRestart = niterrestart; }
  void setEps(double eps) { _eps = eps; }
  void setPrecond(const ALinearOpMulti* precond, Id status);

  const LogStats& getLogStats() const { return _logStats; }

  void prepare() const;

  void setUserInitialValue(bool b) { _userInitialValue = b; }

#ifndef SWIG

protected:
  virtual void _evalDirect(const std::vector<std::vector<double>>& inv,
                           std::vector<std::vector<double>>& outv) const = 0;

public:
  void evalDirect(const std::vector<std::vector<double>>& inv,
                  std::vector<std::vector<double>>& outv) const;
  virtual void evalInverse(const std::vector<std::vector<double>>& vecin,
                           std::vector<std::vector<double>>& vecout) const;
#endif

protected:
  void _updated() const;

private:
  Id                       _nIterMax;
  Id                       _nIterRestart;
  double                    _eps;
  bool                      _precondStatus;
  bool                      _userInitialValue;
  const ALinearOpMulti*     _precond;

  // Work arrays
  mutable bool                         _initialized;
  mutable std::vector<std::vector<double>> _r;

public:
  mutable std::vector<std::vector<double>> _temp;
  mutable std::vector<std::vector<double>> _p;
  mutable std::vector<std::vector<double>> _z;
  mutable double _nb;

protected:
  LogStats                   _logStats;
};
}