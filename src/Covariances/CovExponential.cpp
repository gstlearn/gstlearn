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
