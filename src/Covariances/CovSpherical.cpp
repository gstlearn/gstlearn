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
#include "Covariances/CovSpherical.hpp"

#include "Covariances/CovContext.hpp"

CovSpherical::CovSpherical(const CovContext& ctxt)
: ACovFunc(COV_SPHERICAL, ctxt)
{
}

CovSpherical::CovSpherical(const CovSpherical &r)
: ACovFunc(r)
{
}

CovSpherical& CovSpherical::operator=(const CovSpherical &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovSpherical::~CovSpherical()
{
}

double CovSpherical::_evaluateCov(double h) const
{
  double cov = 0.;
  if (h < 1) cov = 1 - 0.5 * h * (3 - h * h);
  return (cov);
}

String CovSpherical::getFormula() const
{
  return "C(h)=1-\\frac{3}{2}\\left(\\frac{h}{a}\\right)+ \\frac{1}{2}\\left(\\frac{h}{a}\\right)^3";
}
