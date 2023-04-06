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

class GSTLEARN_EXPORT CovReg1D : public ACovFunc
{
public:
  CovReg1D(const CovContext& ctx);
  CovReg1D(const CovReg1D &r);
  CovReg1D& operator= (const CovReg1D &r);
  virtual ~CovReg1D();

  virtual String getFormula() const override { return String("Equation not yet implemented"); }
  String         getCovName() const override { return "1-D Regularized"; }
  int            getMinOrder() const override { return -1; }

  unsigned int getMaxNDim()   const  override { return 1; }
  double getScadef() const override;

protected:
  double _evaluateCov(double h) const override;
};

