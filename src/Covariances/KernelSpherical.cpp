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
#include "Covariances/KernelSpherical.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"

namespace gstlrn
{
KernelSpherical::KernelSpherical(const CovContext& ctxt)
  : AKernel(ECov::SPHERICAL, ctxt)
{
}

KernelSpherical::KernelSpherical(const KernelSpherical& r)
  : AKernel(r)
{
}

KernelSpherical& KernelSpherical::operator=(const KernelSpherical& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelSpherical::~KernelSpherical()
{
}

double KernelSpherical::_evaluateCov(double h) const
{
  double cov = 0.;
  if (h < 1) cov = 1 - 0.5 * h * (3. - h * h);
  return (cov);
}

double KernelSpherical::_evaluateCovFirstDerivative(double h) const
{
  double cov = 0.;
  if (h < 1) cov = -1.5 + 1.5 * h * h;
  return (cov);
}

String KernelSpherical::getFormula() const
{
  return "C(h)=1-\\frac{3}{2}\\left(\\frac{h}{a}\\right)+ \\frac{1}{2}\\left(\\frac{h}{a}\\right)^3";
}

double KernelSpherical::simulateTurningBand(double t0, TurningBandOperate& operTB) const
{
  return operTB.shotNoiseAffineOne(t0);
}
}