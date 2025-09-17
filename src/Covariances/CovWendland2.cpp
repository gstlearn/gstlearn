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
#include "Covariances/CovWendland2.hpp"
#include "Covariances/CovContext.hpp"

namespace gstlrn
{
CovWendland2::CovWendland2(const CovContext& ctxt)
  : ACovFunc(ECov::WENDLAND2, ctxt)
{
}

CovWendland2::CovWendland2(const CovWendland2& r)
  : ACovFunc(r)
{
}

CovWendland2& CovWendland2::operator=(const CovWendland2& r)
{
  if (this != &r)
  {
    ACovFunc::operator=(r);
  }
  return *this;
}

CovWendland2::~CovWendland2()
{
}

double CovWendland2::_evaluateCov(double h) const
{
  // From "Computed Supported Correlation Functions" by T. Gneiting with n=4
  double cov = 0.;
  double h2  = h * h;
  if (h < 1)
    cov = (1 + 6. * h + h2 * (35. / 3.)) * pow(1 - h, 6.);
  return (cov);
}

double CovWendland2::_evaluateCovFirstDerivativeOverH(double h) const
{
  double res;

  res = 0.;
  if (h > 1) return res;
  res = (-(56. / 3.) - h * (280. / 3.)) * pow(1 - h, 5);
  return res;
}

double CovWendland2::_evaluateCovDerivative(Id degree, double h) const
{
  double h2, res;

  res = 0.;
  h2  = h * h;
  if (h > 1) return res;

  switch (degree)
  {
    case 1:
      res = (-h * 56. - h2 * 280.) * pow(1 - h, 5) / 3.;
      break;

    case 2:
      res = (-56. + 784. * h + 840. * h2) * pow(1 - h, 4.);
      break;
  }
  return (res);
}
} // namespace gstlrn