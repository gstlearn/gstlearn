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

class GSTLEARN_EXPORT CovExponential : public ACovFunc
{
public:
  CovExponential(const CovContext& ctx);
  CovExponential(const CovExponential &r);
  CovExponential& operator= (const CovExponential &r);
  virtual ~CovExponential();

  double getScadef() const override;
  virtual String getFormula() const override;
  String         getCovName() const override { return "Exponential"; }
  int            getMinOrder() const override { return -1; }

protected:
  double _evaluateCov(double h) const override;
};

