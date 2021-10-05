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
