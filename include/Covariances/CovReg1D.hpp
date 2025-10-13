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
class CovContext;

class GSTLEARN_EXPORT CovReg1D : public ACovFunc
{
public:
  CovReg1D(const CovContext& ctx);
  CovReg1D(const CovReg1D &r);
  CovReg1D& operator= (const CovReg1D &r);
  virtual ~CovReg1D();

  double         getScadef() const override;
  String         getCovName() const override { return "1-D Regularized"; }
  Id            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

  size_t getMaxNDim() const override { return 1; }

protected:
  double _evaluateCov(double h) const override;
};
}
