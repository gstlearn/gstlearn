/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

class CovContext;

class GSTLEARN_EXPORT CovCauchy : public ACovFunc
{
public:
  CovCauchy(const CovContext& ctx);
  CovCauchy(const CovCauchy &r);
  CovCauchy& operator= (const CovCauchy &r);
  virtual ~CovCauchy();

  virtual String getFormula() const override;
  String         getCovName() const override { return "Cauchy"; }
  int            getMinOrder() const override { return -1; }

  bool   hasParam()  const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;

protected:
  double _evaluateCov(double h) const override;
};

