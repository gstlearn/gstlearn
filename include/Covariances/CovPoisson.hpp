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

class CovContext;

class GSTLEARN_EXPORT CovPoisson : public ACovFunc
{
public:
  CovPoisson(const CovContext& ctx);
  CovPoisson(const CovPoisson &r);
  CovPoisson& operator= (const CovPoisson &r);
  virtual ~CovPoisson();

  String         getCovName() const override { return "Poisson"; }
  bool           hasParam() const override { return true; }
  double         getParMax() const override { return MAX_PARAM; }
  int            getMinOrder() const override { return -1; }
  bool           getCompatibleSpaceS() const override { return true; }
  bool           hasCovOnSphere() const override { return true; }
  bool           hasSpectrumOnSphere() const override { return true; }

protected:
  double _evaluateCovOnSphere(double alpha,
                              double scale = 1.,
                              int degree = 50) const override;
  VectorDouble _evaluateSpectrumOnSphere(int n, double scale = 1.) const override;
};

