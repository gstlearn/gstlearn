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

#include "Covariances/ACovFunc.hpp"
#include "gstlearn_export.hpp"

// In Piecewise polynomial, positive definite and compactly supported
// radial functions of minimal degree, by H. Wendland
// Advances in Computational Mathematics, Vol. 4 (389-396), 1995
// It corresponds to Wendland \psi_{4,2}

namespace gstlrn
{
class CovContext;

class GSTLEARN_EXPORT CovWendland2: public ACovFunc
{
public:
  CovWendland2(const CovContext& ctx);
  CovWendland2(const CovWendland2& r);
  CovWendland2& operator=(const CovWendland2& r);
  virtual ~CovWendland2();

  size_t getMaxNDim() const override { return 3; }
  String getCovName() const override { return "Wendland-4,2"; }
  Id getMinOrder() const override { return -1; }
  bool getCompatibleSpaceR() const override { return true; }
  bool hasCovDerivative() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovDerivative(Id degree, double h) const override;
  double _evaluateCovFirstDerivativeOverH(double h) const override;
};

} // namespace gstlrn