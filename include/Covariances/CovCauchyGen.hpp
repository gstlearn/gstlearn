/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Covariances/ACovFunc.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class CovContext;

class GSTLEARN_EXPORT CovCauchyGen: public ACovFunc
{
public:
  CovCauchyGen(const CovContext& ctx);
  CovCauchyGen(const CovCauchyGen& r);
  CovCauchyGen& operator=(const CovCauchyGen& r);
  virtual ~CovCauchyGen();

  virtual String getFormula() const override;
  String getCovName() const override { return "Generalized Cauchy"; };
  Id getMinOrder() const override { return -1; };
  bool getCompatibleSpaceR() const override { return true; };

  bool hasParam() const override { return true; };
  double getParMax() const override { return MAX_PARAM; };
  double getScadef() const override;
  double getAlpha() const { return _alpha; };
  bool setAlpha(double alpha);

protected:
  double _evaluateCov(double h) const override;

private:
  double _alpha;
};

} // namespace gstlrn