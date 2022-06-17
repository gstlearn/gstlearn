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

class CovContext;

class GSTLEARN_EXPORT CovGaussian : public ACovFunc
{
public:
  CovGaussian(const CovContext& ctx);
  CovGaussian(const CovGaussian &r);
  CovGaussian& operator= (const CovGaussian &r);
  virtual ~CovGaussian();

  virtual String getFormula() const override;
  String         getCovName() const override { return "Gaussian"; }
  int            getMinOrder() const override { return -1; }
  double         getScadef() const override;
  virtual bool   hasCovDerivative() const override { return true; }

protected:
  double _evaluateCov(double h)  const override;
  double _evaluateCovDerivate(int degree, double h) const override;
};

