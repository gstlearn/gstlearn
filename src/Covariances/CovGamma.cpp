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
#include "Covariances/CovGamma.hpp"

#include "math.h"
#include "Covariances/CovContext.hpp"

CovGamma::CovGamma(const CovContext& ctxt)
: ACovFunc(COV_GAMMA, ctxt)
{
  setParam(1);
}

CovGamma::CovGamma(const CovGamma &r)
: ACovFunc(r)
{
}

CovGamma& CovGamma::operator=(const CovGamma &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovGamma::~CovGamma()
{
}

double CovGamma::getScadef() const
{
  return (pow(20.,1. / getParam()) - 1.);
}

double CovGamma::_evaluateCov(double h) const
{
  double cov;
  cov = 1. / pow(1. + h, getParam());
  return (cov);
}

String CovGamma::getFormula() const
{
  return "C(h)=\\frac{1}{\\left( 1+ \\frac{h}{a_t} \\right)^\\alpha";
}
