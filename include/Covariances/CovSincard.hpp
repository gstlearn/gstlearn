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

class CovSincard : public ACovFunc
{
public:
  CovSincard(const CovContext& ctx);
  CovSincard(const CovSincard &r);
  CovSincard& operator= (const CovSincard &r);
  virtual ~CovSincard();

  double getScadef() const override;
  virtual String getFormula() const override;
  String         getCovName() const override { return "Cardinal Sine"; }

protected:
  double _evaluateCov(double h) const override;
};

