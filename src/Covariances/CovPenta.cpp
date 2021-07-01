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
#include "Covariances/CovPenta.hpp"

#include "Covariances/CovContext.hpp"

CovPenta::CovPenta(const CovContext& ctxt)
: ACovFunc(COV_PENTA, ctxt)
{
}

CovPenta::CovPenta(const CovPenta &r)
: ACovFunc(r)
{
}

CovPenta& CovPenta::operator=(const CovPenta &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovPenta::~CovPenta()
{
}

double CovPenta::_evaluateCov(double h) const
{
  double cov = 0.;
  if (h < 1.)
  {
    cov = 1. - 3. * h * (1. - h / 2. * (1. + h / 6.));
  }
  else if (h < 2.)
  {
    cov = -2. + 3. * h * (1. - h / 2. * (1. - h / 6.));
  }
  return (cov);
}

