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

#include "Basic/VectorNumT.hpp"
#include "LinearOp/LogStats.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
class LogStats;

class GSTLEARN_EXPORT ALinearOpMulti
{

public:
  ALinearOpMulti(Id nitermax = 1000, double eps = EPSILON8);
  ALinearOpMulti(const ALinearOpMulti& m);
  ALinearOpMulti& operator=(const ALinearOpMulti& m);
  virtual ~ALinearOpMulti();

  void initLk(const VectorVectorDouble& inv, VectorVectorDouble& outv) const;
  virtual Id sizes() const  = 0;
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
  virtual void _evalDirect(const VectorVectorDouble& inv,
                           VectorVectorDouble& outv) const = 0;

public:
  void evalDirect(const VectorVectorDouble& inv,
                  VectorVectorDouble& outv) const;
  virtual void evalInverse(const VectorVectorDouble& vecin,
                           VectorVectorDouble& vecout) const;
#endif

protected:
  void _updated() const;

private:
  Id _nIterMax;
  Id _nIterRestart;
  double _eps;
  bool _precondStatus;
  bool _userInitialValue;
  const ALinearOpMulti* _precond;

  // Work arrays
  mutable bool _initialized;
  mutable VectorVectorDouble _r;

public:
  mutable VectorVectorDouble _temp;
  mutable VectorVectorDouble _p;
  mutable VectorVectorDouble _z;
  mutable double _nb;

protected:
  LogStats _logStats;
};
} // namespace gstlrn