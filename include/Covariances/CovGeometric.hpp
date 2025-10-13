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

class GSTLEARN_EXPORT CovGeometric : public ACovFunc
{
public:
  CovGeometric(const CovContext& ctx);
  CovGeometric(const CovGeometric &r);
  CovGeometric& operator= (const CovGeometric &r);
  virtual ~CovGeometric();

  String         getCovName() const override { return "Geometric"; }
  Id            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceS() const override { return true; }
  bool           hasCovOnSphere() const override { return true; }
  bool           hasSpectrumOnSphere() const override { return true; }

protected:
  double _evaluateCovOnSphere(double alpha,
                              double scale = 1.,
                              Id degree = 50) const override;
  VectorDouble _evaluateSpectrumOnSphere(Id n, double scale = 1.) const override;
};

}