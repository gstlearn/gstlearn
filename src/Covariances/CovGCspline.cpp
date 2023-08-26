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
#include "Covariances/CovGCspline.hpp"

#include <math.h>
#include "Covariances/CovContext.hpp"

CovGCspline::CovGCspline(const CovContext& ctxt)
: ACovFunc(ECov::SPLINE_GC, ctxt)
{
}

CovGCspline::CovGCspline(const CovGCspline &r)
: ACovFunc(r)
{
}

CovGCspline& CovGCspline::operator=(const CovGCspline &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovGCspline::~CovGCspline()
{
}

double CovGCspline::_evaluateCov(double h) const
{
  int ndim = getContext().getNDim();
  double r = getContext().getField();
  double r2 = r * r;
  double h2 = h * h;
  double logval = (r < 10.e-5 || h < r * 10.e-5) ? 0. : log(h / r);

  // Code for the first 3 Space dimensions
  double cov = 0.;
  if (ndim == 1)
    cov = 0.5 * r2 - h2 * (1.5 - log(2) - logval);
  else if (ndim == 2)
    cov = r2 - h2 * (1 - logval);
  else
    cov = 1.5 * r2 - h2 * (11. / 6. - log(2) - logval);

  return (cov);
}
