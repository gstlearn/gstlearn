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
#include "Covariances/KernelGamma.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
KernelGamma::KernelGamma(const CovContext& ctxt)
  : AKernel(ECov::GAMMA, ctxt)
{
  setParam(1);
}

KernelGamma::KernelGamma(const KernelGamma& r)
  : AKernel(r)
{
}

KernelGamma& KernelGamma::operator=(const KernelGamma& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelGamma::~KernelGamma()
{
}

double KernelGamma::getScadef() const
{
  double param = getParam();
  if (param < 0.05) param = 1.; // This test is too avoid passing too small number to next line
  double scadef = pow(20., 1. / param) - 1.;
  return (scadef);
}

double KernelGamma::_evaluateCov(double h) const
{
  double cov;
  cov = 1. / pow(1. + h, getParam());
  return (cov);
}

String KernelGamma::getFormula() const
{
  return "C(h)=\\frac{1}{\\left( 1+ \\frac{h}{a_t} \\right)^\\alpha";
}
} // namespace gstlrn