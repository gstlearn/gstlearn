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

class ALinearOp;

class GSTLEARN_EXPORT CGParam {

public:
  CGParam(int nitermax = 1000, double eps = EPSILON8);
  CGParam(const CGParam &m);
  CGParam& operator=(const CGParam &m);
  virtual ~CGParam();

  void setEps(double eps) { _eps = eps; }
  void setNIterMax(int nIterMax) { _nIterMax = nIterMax; }
  void setX0(VectorDouble x0) { _x0 = x0; }
  void setPrecond(const ALinearOp* precond, int status);
  void setPrecondStatus(int precondStatus) { _precondStatus = precondStatus; }

  double getEps() const { return _eps; }
  int getNIterMax() const { return _nIterMax; }
  const ALinearOp* getPrecond() const { return _precond; }
  const VectorDouble& getX0() const { return _x0; }
  double getX0(int i) const { return _x0[i]; }
  int getPrecondStatus() const { return _precondStatus; }

private:
  int              _nIterMax;
  double           _eps;
  VectorDouble     _x0;
  int              _precondStatus;
  const ALinearOp* _precond; // External Pointer
};
