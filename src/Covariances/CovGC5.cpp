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
#include "Covariances/CovGC5.hpp"

#include "Covariances/CovContext.hpp"

CovGC5::CovGC5(const CovContext& ctxt)
: ACovFunc(ECov::ORDER5_GC, ctxt)
{
}

CovGC5::CovGC5(const CovGC5 &r)
: ACovFunc(r)
{
}

CovGC5& CovGC5::operator=(const CovGC5 &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovGC5::~CovGC5()
{
}

double CovGC5::_evaluateCov(double h) const
{
  double cov;
  double r = getContext().getField();
  int ndim = getContext().getNDim();
  double h2 = h * h;
  double r2 = r * r;
  double h4 = h2 * h2;
  double r3 = r2 * r;

  if (ndim == 1)
    cov = h4 * (h - 5. * r) + r3 * (20. * h2 - 16. * r2);
  else if (ndim == 2)
    cov = h4 * (h - 225. * GV_PI * r / 128.)
        + r3 * (75. * GV_PI * h2 / 8. - 15. * GV_PI * r2);
  else
    cov = h4 * (h - 6. * r) + r3 * (40. * h2 - 96. * r2);

  return (-cov);
}
