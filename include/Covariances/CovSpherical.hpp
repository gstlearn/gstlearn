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

class CovSpherical : public ACovFunc
{
public:
  CovSpherical(const CovContext& ctx);
  CovSpherical(const CovSpherical &r);
  CovSpherical& operator= (const CovSpherical &r);
  virtual ~CovSpherical();

  unsigned int getMaxNDim() const override { return 3; }
  String       getFormula() const override;
  String       getCovName() const override { return "Spherical"; }

protected:
  double _evaluateCov(double h) const override;
};

