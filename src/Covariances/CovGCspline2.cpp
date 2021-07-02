/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Covariances/CovGCspline2.hpp"

#include <math.h>
#include "Covariances/CovContext.hpp"

CovGCspline2::CovGCspline2(const CovContext& ctxt)
: ACovFunc(COV_SPLINE2_GC, ctxt)
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

double CovGCspline2::_evaluateCovDerivate(int degree, double h) const
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
