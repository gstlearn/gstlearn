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

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"


namespace gstlrn
{
class ALinearOp;

class GSTLEARN_EXPORT CGParam: public AStringable
{
public:
  CGParam(Id nitermax = 1000, double eps = EPSILON8);
  CGParam(const CGParam& m);
  CGParam& operator=(const CGParam& m);
  virtual ~CGParam();

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  void setEps(double eps) { _eps = eps; }
  void setNIterMax(Id nIterMax) { _nIterMax = nIterMax; }
  void setX0(const VectorDouble& x0) { _x0 = x0; }
  void setPrecond(const ALinearOp* precond, Id status);
  void setPrecondStatus(Id precondStatus) { _precondStatus = precondStatus; }

  double getEps() const { return _eps; }
  Id getNIterMax() const { return _nIterMax; }
  Id getPrecondStatus() const { return _precondStatus; }
  double getX0(Id i) const { return _x0[i]; }
  const ALinearOp* getPrecond() const { return _precond; }
  const VectorDouble& getX0() const { return _x0; }

private:
  Id _nIterMax;
  double _eps;
  VectorDouble _x0;
  Id _precondStatus;
  const ALinearOp* _precond; // External Pointer
};
}