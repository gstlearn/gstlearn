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

class GSTLEARN_EXPORT CovCosExp : public ACovFunc
{
public:
  CovCosExp(const CovContext& ctx);
  CovCosExp(const CovCosExp &r);
  CovCosExp& operator= (const CovCosExp &r);
  virtual ~CovCosExp();

  double getParMax() const override { return TEST; }
  bool   hasParam()  const override { return true; }
  double getScadef() const override;

  String getCovName() const override { return "Cosexp"; }

protected:
  double _evaluateCov(double h) const override;
};

