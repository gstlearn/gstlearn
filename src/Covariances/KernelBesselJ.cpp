/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Covariances/KernelBesselJ.hpp"
#include "Basic/MathFunc.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"

#include <cmath>

#define MAXTAB 100

namespace gstlrn
{

KernelBesselJ::KernelBesselJ(const CovContext& ctxt)
  : AKernel(ECov::BESSELJ, ctxt)
{
  setParam(1);
}

KernelBesselJ::KernelBesselJ(const KernelBesselJ& r)
  : AKernel(r)
{
}

KernelBesselJ& KernelBesselJ::operator=(const KernelBesselJ& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelBesselJ::~KernelBesselJ()
{
}

double KernelBesselJ::_evaluateCov(double h) const
{
  static double TAB[MAXTAB];

  double cov   = 0.;
  double third = getParam();
  Id nb        = static_cast<Id>(floor(third));
  double alpha = third - nb;
  if (third <= 0 || nb >= MAXTAB) return (cov);
  double coeff = (h > 0) ? pow(h / 2., third) : 1.;

  cov = 1.;
  if (h > 0)
  {
    if (besselj_table(h, alpha, nb + 1, TAB) < nb + 1) return (cov);
    cov = TAB[nb] * exp(loggamma(third + 1.)) / coeff;
  }
  return (cov);
}

String KernelBesselJ::getFormula() const
{
  return "C(h)=2^\\alpha\\Gamma(\\alpha+1) \\frac{ J_\\alpha\\left( \\frac{h}{a_t} \\right) } {\\left( \\frac{h}{a_t} \\right)^\\alpha}";
}

double KernelBesselJ::simulateTurningBand(double t0, TurningBandOperate& operTB) const
{
  return operTB.cosineOne(t0);
}
} // namespace gstlrn