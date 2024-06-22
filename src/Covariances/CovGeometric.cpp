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
#include "Covariances/CovGeometric.hpp"

#include "math.h"
#include "Basic/VectorHelper.hpp"
#include "Covariances/CovContext.hpp"

CovGeometric::CovGeometric(const CovContext& ctxt)
: ACovFunc(ECov::GEOMETRIC, ctxt)
{
  setParam(1);
}

CovGeometric::CovGeometric(const CovGeometric &r)
: ACovFunc(r)
{
}

CovGeometric& CovGeometric::operator=(const CovGeometric &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovGeometric::~CovGeometric()
{
}

double CovGeometric::_evaluateCovOnSphere(double alpha,
                                          double scale,
                                          double param,
                                          int degree) const
{
  DECLARE_UNUSED(param);
  double rho = scale;
  return ((1. - rho) / sqrt(1. - 2. * rho * cos(alpha) + rho * rho));
}

VectorDouble CovGeometric::_evaluateSpectrumOnSphere(int n,
                                                     double scale,
                                                     double param) const
{
  DECLARE_UNUSED(param);
  double rho = scale;
  VectorDouble sp(1+n, 0.);

  double rhoprod = 1.;
  for (int k = 0; k <= n; k++)
  {
    sp[k] = rhoprod;
    rhoprod *= rho;
  }

  VH::normalize(sp, 1);

  return sp;
}
