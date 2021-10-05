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
#include "Covariances/CovTriangle.hpp"

#include "Covariances/CovContext.hpp"

CovTriangle::CovTriangle(const CovContext& ctxt)
: ACovFunc(ECov::TRIANGLE, ctxt)
{
}

CovTriangle::CovTriangle(const CovTriangle &r)
: ACovFunc(r)
{
}

CovTriangle& CovTriangle::operator=(const CovTriangle &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovTriangle::~CovTriangle()
{
}

double CovTriangle::_evaluateCov(double h) const
{
  double cov = MAX(0, 1. - h);
  return (cov);
}

