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

class GSTLEARN_EXPORT CovP8 : public ACovFunc
{
public:
  CovP8(const CovContext& ctx);
  CovP8(const CovP8 &r);
  CovP8& operator= (const CovP8 &r);
  virtual ~CovP8();

  unsigned int getMaxNDim()   const  override { return 3; }

  virtual String getFormula() const override { return String("Equation not yet implemented"); }
  String         getCovName() const override { return "Order-8 Polynomial"; }
  virtual bool   hasCovDerivative() const override { return true; }

protected:
  double _evaluateCov(double h)  const override;
  double _evaluateCovDerivate(int degree, double h) const override;
};

