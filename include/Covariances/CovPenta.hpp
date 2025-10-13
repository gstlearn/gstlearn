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

class GSTLEARN_EXPORT CovPenta : public ACovFunc
{
public:
  CovPenta(const CovContext& ctx);
  CovPenta(const CovPenta &r);
  CovPenta& operator= (const CovPenta &r);
  virtual ~CovPenta();

  size_t getMaxNDim() const override { return 3; }

  String         getCovName() const override { return "Penta"; }
  Id            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceR() const override { return true; }

protected:
  double _evaluateCov(double h) const override;
};
}
