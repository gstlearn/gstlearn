/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

class CovContext;

class GSTLEARN_EXPORT CovGamma : public ACovFunc
{
public:
  CovGamma(const CovContext& ctx);
  CovGamma(const CovGamma &r);
  CovGamma& operator= (const CovGamma &r);
  virtual ~CovGamma();

  virtual String getFormula() const override;
  String         getCovName() const override { return "Gamma"; }
  int            getMinOrder() const override { return -1; }

  bool   hasParam() const override { return true; }
  double getParMax() const override { return MAX_PARAM; }
  double getScadef() const override;

protected:
  double _evaluateCov(double h) const override;
};

