/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Covariances/CovGamma.hpp"

#include "math.h"
#include "Covariances/CovContext.hpp"

CovGamma::CovGamma(const CovContext& ctxt)
: ACovFunc(ECov::GAMMA, ctxt)
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
