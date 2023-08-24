/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

/* Be careful ! This is not a real covariance */

class CovContext;

class GSTLEARN_EXPORT CovGCspline2 : public ACovFunc
{
public:
  CovGCspline2(const CovContext& ctx);
  CovGCspline2(const CovGCspline2 &r);
  CovGCspline2& operator= (const CovGCspline2 &r);
  virtual ~CovGCspline2();

  int          getMinOrder()  const override { return 1; }
  unsigned int getMaxNDim()   const  override { return 3; }
  String       getCovName() const override { return "Spline-2 G.C."; }
  virtual bool hasCovDerivative() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovDerivative(int degree, double h) const override;
};

