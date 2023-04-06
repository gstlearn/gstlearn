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
#include "Covariances/CovPower.hpp"
#include "Covariances/CovContext.hpp"
#include "Basic/MathFunc.hpp"

#include <math.h>

CovPower::CovPower(const CovContext& ctxt)
: ACovFunc(ECov::POWER, ctxt)
{
  setParam(1);
}

CovPower::CovPower(const CovPower &r)
: ACovFunc(r)
{
}

CovPower& CovPower::operator=(const CovPower &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovPower::~CovPower()
{
}

double CovPower::_evaluateCov(double h) const
{
  double a;

  int ndim = getContext().getNDim();
  double alpha = getParam();
  double r = getContext().getField();
  double ra = pow(r, alpha);
  double g1 = loggamma((alpha + 1.) / 2.) + loggamma(1. - alpha / 2.);

  if (ndim == 1)
    a = exp(g1) * ra / sqrt(GV_PI);
  else if (ndim == 2)
    a = (exp(loggamma(alpha + 1.5) - loggamma(alpha + 1.) + g1) * 2. * ra / GV_PI);
  else
    a = (alpha + 1.) * exp(g1) * ra / sqrt(GV_PI);

  double cov = a;
  if (h > 0) cov -= pow(h, getParam());
  return (cov);
}

