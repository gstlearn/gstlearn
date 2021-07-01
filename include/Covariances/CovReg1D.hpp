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

#include "Covariances/ACovFunc.hpp"

class CovContext;

class CovReg1D : public ACovFunc
{
public:
  CovReg1D(const CovContext& ctx);
  CovReg1D(const CovReg1D &r);
  CovReg1D& operator= (const CovReg1D &r);
  virtual ~CovReg1D();

  virtual String getFormula() const override { return String("Equation not yet implemented"); }
  String         getCovName() const override { return "1-D Regularized"; }

  unsigned int getMaxNDim()   const  override { return 1; }
  double getScadef() const override;

protected:
  double _evaluateCov(double h) const override;
};

