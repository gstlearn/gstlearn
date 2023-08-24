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
#include "Covariances/CovGCspline2.hpp"

#include <math.h>
#include "Covariances/CovContext.hpp"

CovGCspline2::CovGCspline2(const CovContext& ctxt)
: ACovFunc(ECov::SPLINE2_GC, ctxt)
{
}

CovGCspline2::CovGCspline2(const CovGCspline2 &r)
: ACovFunc(r)
{
}

CovGCspline2& CovGCspline2::operator=(const CovGCspline2 &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovGCspline2::~CovGCspline2()
{
}

double CovGCspline2::_evaluateCov(double h) const
{
  double B = 1.;
  double A = (7. - 10. * B) / 12.;
  double C = (-7. - 2. * B) / 12.;

  double h2 = h * h;
  double logval = (h < 10.e-5) ? 0. : log(h);
  double cov = A + h2 * (B + h2 * (C + logval));

  return (cov);
}

double CovGCspline2::_evaluateCovDerivative(int degree, double h) const
{
  double B = 1.;
  double C = (-7. - 2. * B) / 12.;
  double h2 = h * h;
  double logval = (h < 10.e-5) ? 0. : log(h);

  double cov = 0.;
  switch (degree)
  {
    case 1:
      cov = 2. * B + h2 * (1. + 4. * (C + logval));
      break;

    case 2:
      cov = 2. * B + h2 * (7. + 12. * (C + logval));
      break;
  }

  return (cov);
}
