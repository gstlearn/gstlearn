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
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT ALinearOp {

public:
  ALinearOp(int nitermax = 1000, double eps = EPSILON8);
  ALinearOp(const ALinearOp &m);
  ALinearOp& operator=(const ALinearOp &m);
  virtual ~ALinearOp();

  virtual void evalInverse(const VectorDouble& inv, VectorDouble& outv) const;
  virtual int getSize() const = 0;

  void evalDirect(const VectorDouble& inv, VectorDouble& outv) const;
  void setNIterMax(int nitermax) { _nIterMax = nitermax; }
  void setEps(double eps) { _eps = eps; }
  void setX0(VectorDouble& x0) { _x0 = x0; }
  void setPrecond(const ALinearOp* precond, int status);

protected:
  virtual void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const = 0;

private:
  double _prod(const VectorDouble& x, const VectorDouble& y) const;

private:
  int              _nIterMax;
  double           _eps;
  VectorDouble     _x0;
  int              _precondStatus;
  const ALinearOp* _precond; // Pointer copied
};
