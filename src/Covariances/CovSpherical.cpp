/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovSpherical.hpp"

#include "Covariances/CovContext.hpp"

CovSpherical::CovSpherical(const CovContext& ctxt)
: ACovFunc(ECov::SPHERICAL, ctxt)
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
