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
#include "Covariances/CovBesselK.hpp"
#include "Covariances/CovContext.hpp"
#include "Basic/MathFunc.hpp"
#include "geoslib_f.h"
#include "math.h"

#define MAXTAB 100

CovBesselK::CovBesselK(const CovContext& ctxt)
: ACovFunc(ECov::BESSEL_K, ctxt)
{
  setParam(1);
}

CovBesselK::CovBesselK(const CovBesselK &r)
: ACovFunc(r)
{
}

CovBesselK& CovBesselK::operator=(const CovBesselK &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovBesselK::~CovBesselK()
{
}

double CovBesselK::getScadef() const
{
  return sqrt(12. * getParam());
}

double CovBesselK::_evaluateCov(double h) const
{
  static double TAB[MAXTAB];

  double cov = 0.;
  double third = getParam();
  int nb = (int) floor(third);
  double alpha = third - nb;
  if (third <= 0 || nb >= MAXTAB) return (0.);
  double coeff = (h > 0) ? pow(h / 2., third) :
                    1.;
  cov = 1.;
  if (h > 0)
  {
    if (bessel_k(h, alpha, nb + 1, TAB) < nb + 1) return (cov);
    cov = 2. * coeff * TAB[nb] / exp(loggamma(third));
  }
  return (cov);
}

String CovBesselK::getFormula() const
{
  return "C(h)=\\frac{ \\left( \\frac{h}{a_t} \\right)^\\alpha}{2^{\\alpha-1}\\Gamma(\\alpha)}K_{-\\alpha} \\left( \\frac{h}{a_t} \\right)";
}

double CovBesselK::_evaluateCovOnSphere(double scale, int degree) const
{
  double kappa = 1. / scale;
  return (2.*degree+1.) / pow(kappa * kappa + degree * (degree+1), getParam());
}
