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

class CovCubic : public ACovFunc
{
public:
  CovCubic(const CovContext& ctx);
  CovCubic(const CovCubic &r);
  CovCubic& operator= (const CovCubic &r);
  virtual ~CovCubic();

  unsigned int getMaxNDim()   const  override { return 3; }

  virtual String getFormula() const override;
  String         getCovName() const override { return "Cubic"; }

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovDerivate(int degree, double h) const override;
};

