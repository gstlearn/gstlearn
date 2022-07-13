/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

// In Piecewise polynomial, positive definite and compactly supported
// radial functions of minimal degree, by H. Wendland
// Advances in Computational Mathematics, Vol. 4 (389-396), 1995
// It corresponds to Wendland \psi_{2,0}

class CovContext;

class GSTLEARN_EXPORT CovWendland0 : public ACovFunc
{
public:
  CovWendland0(const CovContext& ctx);
  CovWendland0(const CovWendland0 &r);
  CovWendland0& operator= (const CovWendland0 &r);
  virtual ~CovWendland0();

  unsigned int getMaxNDim()   const  override { return 3; }

  virtual String getFormula() const override { return String("Equation not yet implemented"); }
  String         getCovName() const override { return "Wendland-2,0"; }
  int            getMinOrder() const override { return -1; }

protected:
  double _evaluateCov(double h) const override;
};

