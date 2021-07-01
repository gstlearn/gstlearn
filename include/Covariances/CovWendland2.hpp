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

class CovWendland2 : public ACovFunc
{
public:
  CovWendland2(const CovContext& ctx);
  CovWendland2(const CovWendland2 &r);
  CovWendland2& operator= (const CovWendland2 &r);
  virtual ~CovWendland2();

  unsigned int getMaxNDim()   const  override { return 3; }

  virtual String getFormula() const override { return String("Equation not yet implemented"); }
  String         getCovName() const override { return "Wendland-2"; }

protected:
  double _evaluateCov(double h) const override;
};

