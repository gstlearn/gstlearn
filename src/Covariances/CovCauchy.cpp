/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovCauchy.hpp"
#include "Covariances/CovContext.hpp"

#include "math.h"

CovCauchy::CovCauchy(const CovContext& ctxt)
: ACovFunc(ECov::CAUCHY, ctxt)
{
  setParam(1);
}

CovCauchy::CovCauchy(const CovCauchy &r)
: ACovFunc(r)
{
}

CovCauchy& CovCauchy::operator=(const CovCauchy &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovCauchy::~CovCauchy()
{
}

double CovCauchy::getScadef() const
{
  return sqrt(pow(20., 1. / getParam()) - 1.);
}

double CovCauchy::_evaluateCov(double h) const
{
  double cov = 1. / pow(1. + h * h, getParam());
  return (cov);
}

String CovCauchy::getFormula() const
{
  return "C(h)=\\frac{1}{\\left( 1+ \\frac{h^2}{a_t^2} \\right)^\\alpha";
}
