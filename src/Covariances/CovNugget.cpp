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
#include "Covariances/CovNugget.hpp"

#include "Covariances/CovContext.hpp"

CovNugget::CovNugget(const CovContext& ctxt)
: ACovFunc(COV_NUGGET, ctxt)
{
}

CovNugget::CovNugget(const CovNugget &r)
: ACovFunc(r)
{
}

CovNugget& CovNugget::operator=(const CovNugget &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovNugget::~CovNugget()
{
}

double CovNugget::_evaluateCov(double h) const
{
  double cov = 0.;
  if (ABS(h) < 1.e-10) cov = 1.;
  return (cov);
}

String CovNugget::getFormula() const
{
  return "C(h)=\\delta(h)";
}
