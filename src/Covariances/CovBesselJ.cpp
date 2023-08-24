/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovBesselJ.hpp"
#include "Covariances/CovContext.hpp"
#include "Basic/MathFunc.hpp"

#include <math.h>

#define MAXTAB 100

CovBesselJ::CovBesselJ(const CovContext& ctxt)
: ACovFunc(ECov::BESSEL_J, ctxt)
{
  setParam(1);
}

CovBesselJ::CovBesselJ(const CovBesselJ &r)
: ACovFunc(r)
{
}

CovBesselJ& CovBesselJ::operator=(const CovBesselJ &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovBesselJ::~CovBesselJ()
{
}

double CovBesselJ::_evaluateCov(double h) const
{
  static double TAB[MAXTAB];

  double cov = 0.;
  double third = getParam();
  int nb = (int) floor(third);
  double alpha = third - nb;
  if (third <= 0 || nb >= MAXTAB) return (cov);
  double coeff = (h > 0) ? pow(h / 2., third) : 1.;

  cov = 1.;
  if (h > 0)
  {
    if (bessel_j(h, alpha, nb + 1, TAB) < nb + 1) return (cov);
    cov = TAB[nb] * exp(loggamma(third + 1.)) / coeff;
  }
  return (cov);
}

String CovBesselJ::getFormula() const
{
  return "C(h)=2^\\alpha\\Gamma(\\alpha+1) \\frac{ J_\\alpha\\left( \\frac{h}{a_t} \\right) } {\\left( \\frac{h}{a_t} \\right)^\\alpha}";
}
