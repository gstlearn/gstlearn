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

class CovGamma : public ACovFunc
{
public:
  CovGamma(const CovContext& ctx);
  CovGamma(const CovGamma &r);
  CovGamma& operator= (const CovGamma &r);
  virtual ~CovGamma();

  virtual String getFormula() const override;
  String         getCovName() const override { return "Gamma"; }

  bool   hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;

protected:
  double _evaluateCov(double h) const override;
};

