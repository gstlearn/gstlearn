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

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

namespace gstlrn
{
  // Be careful ! This is not a real covariance
  // It is used to compute the gradient of the covariance function
  // with respect to the parameter of the covariance function
  // It is used in the context of optimization of the covariance function
  // parameters.
class CovContext;

class GSTLEARN_EXPORT CovGamma : public ACovFunc
{
public:
  CovGamma(const CovContext& ctx);
  CovGamma(const CovGamma &r);
  CovGamma& operator= (const CovGamma &r);
  virtual ~CovGamma();

  String getFormula() const override;
  String         getCovName() const override { return "Gamma"; }
  Id            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

  bool   hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;

protected:
  double _evaluateCov(double h) const override;
};

}