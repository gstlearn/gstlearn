/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
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

  void evalDirect(const VectorDouble& inv, VectorDouble& outv) const;
  virtual void evalInverse(const VectorDouble& inv, VectorDouble& outv) const;
  virtual int getSize() const = 0;
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
