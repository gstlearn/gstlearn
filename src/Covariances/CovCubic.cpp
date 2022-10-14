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
#include "Covariances/CovCubic.hpp"

#include "Covariances/CovContext.hpp"

CovCubic::CovCubic(const CovContext& ctxt)
: ACovFunc(ECov::CUBIC, ctxt)
{
}

CovCubic::CovCubic(const CovCubic &r)
: ACovFunc(r)
{
}

CovCubic& CovCubic::operator=(const CovCubic &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovCubic::~CovCubic()
{
}

double CovCubic::_evaluateCov(double h) const
{
  double cov = 0.;
  double h2 = h * h;
  if (h < 1) cov = 1. - h2 * (7. + h * (-8.75 + h2 * (3.5 - 0.75 * h2)));
  cov = MAX(0., cov);
  return (cov);
}

double CovCubic::_evaluateCovDerivative(int degree, double h) const
{
  double h2, res;

  res = 0.;
  h2 = h * h;
  if (h2 >= 1) return res;

  switch (degree)
  {
    case 1:
      res = -14. + h * (26.25 - h2 * (17.5 - 5.25 * h2));
      break;

    case 2:
      res = -14. + h * (52.5 + h2 * (-70. + 31.5 * h2)); // TODO a verifier
      break;
  }
  return (res);
}

String CovCubic::getFormula() const
{
  return "C(h)=1 - h^2 * \\left(7 + h * \\left(-8.75 + h^2 * \\left(3.5 - 0.75 * h^2 \\right) \\right) \\right)";
}
