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
#include "Covariances/CovLinearSph.hpp"

#include "math.h"
#include "Basic/VectorHelper.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/Law.hpp"
#include "Covariances/CovContext.hpp"

CovLinearSph::CovLinearSph(const CovContext &ctxt)
    : ACovFunc(ECov::LINEARSPH, ctxt)
{
  setParam(1);
}

CovLinearSph::CovLinearSph(const CovLinearSph &r)
: ACovFunc(r)
{
}

CovLinearSph& CovLinearSph::operator=(const CovLinearSph &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovLinearSph::~CovLinearSph()
{
}

double CovLinearSph::_evaluateCovOnSphere(double alpha,
                                          double scale,
                                          int degree) const
{
  DECLARE_UNUSED(scale);
  DECLARE_UNUSED(degree);
  return 1. - 2. * alpha / GV_PI;
}

VectorDouble CovLinearSph::_evaluateSpectrumOnSphere(int n, double scale) const
{
  DECLARE_UNUSED(scale);
  VectorDouble sp(n + 1, 0.);

  int k = 1;
  sp[k] = 3. / 4.;
  while (1)
  {
    k += 2;
    if (k >= n + 1) break;
    double v = (k - 2.) / (k + 1.);
    sp[k] = (2. * k + 1.) / (2. * k - 3.) * v * v * sp[k - 2];
  }

  VH::normalize(sp, 1);

  return sp;
}
