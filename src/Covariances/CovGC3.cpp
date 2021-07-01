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
#include "Covariances/CovGC3.hpp"

#include "Covariances/CovContext.hpp"

CovGC3::CovGC3(const CovContext& ctxt)
: ACovFunc(COV_ORDER3_GC, ctxt)
{
}

CovGC3::CovGC3(const CovGC3 &r)
: ACovFunc(r)
{
}

CovGC3& CovGC3::operator=(const CovGC3 &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovGC3::~CovGC3()
{
}

double CovGC3::_evaluateCov(double h) const
{
  double cov;
  double r = getContext().getField();
  int ndim = getContext().getNDim();

  double h2 = h * h;
  double r3 = r * r * r;

  if (ndim == 1)
    cov = h2 * (h - 3. * r) + 2. * r3;
  else if (ndim == 2)
    cov = h2 * (h - 9. * GV_PI * r / 8. * r) + 3. * GV_PI * r3 / 2;
  else
    cov = h2 * (h - 4. * r) + 8. * r3;

  return (cov);
}
