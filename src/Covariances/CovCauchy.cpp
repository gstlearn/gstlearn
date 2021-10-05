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
#include "Covariances/CovCauchy.hpp"

#include "math.h"

#include "Covariances/CovContext.hpp"
#include "geoslib_f.h"

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
