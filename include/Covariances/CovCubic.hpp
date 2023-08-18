/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Covariances/ACovFunc.hpp"

class CovContext;

class GSTLEARN_EXPORT CovCubic : public ACovFunc
{
public:
  CovCubic(const CovContext& ctx);
  CovCubic(const CovCubic &r);
  CovCubic& operator= (const CovCubic &r);
  virtual ~CovCubic();

  unsigned int getMaxNDim()   const  override { return 3; }

  virtual String getFormula() const override;
  String         getCovName() const override { return "Cubic"; }
  int            getMinOrder() const override { return -1; }
  virtual bool   hasCovDerivative() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
  double _evaluateCovDerivative(int degree, double h) const override;
};

