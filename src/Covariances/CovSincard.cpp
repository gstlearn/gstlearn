/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Covariances/CovSincard.hpp"

#include "math.h"
#include "Covariances/CovContext.hpp"

CovSincard::CovSincard(const CovContext& ctxt)
: ACovFunc(ECov::SINCARD, ctxt)
{
}

CovSincard::CovSincard(const CovSincard &r)
: ACovFunc(r)
{
}

CovSincard& CovSincard::operator=(const CovSincard &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovSincard::~CovSincard()
{
}

double CovSincard::getScadef() const
{
  return (20.371);
}

double CovSincard::_evaluateCov(double h) const
{
  static double MIN_SIN = 1.e-5;
  double cov = 1.;
  if (h > MIN_SIN) cov = sin(h) / h;
  return (cov);
}

String CovSincard::getFormula() const
{
  return "C(h)=\\frac{sin(\\frac{h}{a})}{\\frac{h}{a}}";
}
