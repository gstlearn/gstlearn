/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

// In Piecewise polynomial, positive definite and compactly supported
// radial functions of minimal degree, by H. Wendland
// Advances in Computational Mathematics, Vol. 4 (389-396), 1995
// It corresponds to Wendland \psi_{3,1}

class CovContext;

class GSTLEARN_EXPORT CovWendland1 : public ACovFunc
{
public:
  CovWendland1(const CovContext& ctx);
  CovWendland1(const CovWendland1 &r);
  CovWendland1& operator= (const CovWendland1 &r);
  virtual ~CovWendland1();

  unsigned int getMaxNDim()   const  override { return 3; }

  virtual String getFormula() const override { return String("Equation not yet implemented"); }
  String         getCovName() const override { return "Wendland-3,1"; }
  int            getMinOrder() const override { return -1; }

protected:
  double _evaluateCov(double h) const override;
};

