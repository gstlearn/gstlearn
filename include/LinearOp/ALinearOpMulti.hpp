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
#include "Basic/VectorNumT.hpp"
#include "LinearOp/LogStats.hpp"

class ALinearOpMulti;

class GSTLEARN_EXPORT ALinearOpMulti {

public:
  ALinearOpMulti(int nitermax = 1000, double eps = EPSILON8);
  ALinearOpMulti(const ALinearOpMulti &m);
  ALinearOpMulti& operator=(const ALinearOpMulti &m);
  virtual ~ALinearOpMulti();

  virtual void evalInverse(const VectorVectorDouble &vecin,
                           VectorVectorDouble &vecout) const;

  void evalDirect(const VectorVectorDouble &inv,
                  VectorVectorDouble &outv) const;
  void initLk(const VectorVectorDouble &inv, VectorVectorDouble &outv) const;
  virtual int sizes() const = 0;
  virtual int size(int) const = 0;

  void setNIterMax(int nitermax) { _nIterMax = nitermax; }
  void setNIterRestart(int niterrestart) { _nIterRestart = niterrestart; }
  void setEps(double eps) { _eps = eps; }
  void setPrecond(const ALinearOpMulti* precond, int status);

  const LogStats& getLogStats() const { return _logStats; }

  void prepare() const;

  void setUserInitialValue(bool b) { _userInitialValue = b; }

protected:
  virtual void _evalDirect(const VectorVectorDouble &inv,
                           VectorVectorDouble &outv) const = 0;
  void _updated() const;

private:
  int                       _nIterMax;
  int                       _nIterRestart;
  double                    _eps;
  bool                      _precondStatus;
  bool                      _userInitialValue;
  const ALinearOpMulti*     _precond;

  // Work arrays
  mutable bool               _initialized;
  mutable VectorVectorDouble _r;

public:
  mutable VectorVectorDouble _temp;
  mutable VectorVectorDouble _p;
  mutable VectorVectorDouble _z;
  mutable double _nb;

protected:
  LogStats                   _logStats;
};
