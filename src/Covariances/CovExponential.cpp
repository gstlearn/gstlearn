/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Covariances/CovExponential.hpp"

#include "math.h"
#include "Covariances/CovContext.hpp"

CovExponential::CovExponential(const CovContext& ctxt)
: ACovFunc(ECov::EXPONENTIAL, ctxt)
{
}

CovExponential::CovExponential(const CovExponential &r)
: ACovFunc(r)
{
}

CovExponential& CovExponential::operator=(const CovExponential &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovExponential::~CovExponential()
{
}

double CovExponential::getScadef() const
{
  return (2.995732);
}

double CovExponential::_evaluateCov(double h) const
{
  if (h > 100) return (0.);
  double cov = exp(-h);
  return (cov);
}

String CovExponential::getFormula() const
{
  return "C(h)=exp \\left( -\\frac{h}{a_t} \\right)";
}
