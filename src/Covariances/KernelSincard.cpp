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
#include "Covariances/KernelSincard.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"

#include <cmath>

namespace gstlrn
{
  KernelSincard::KernelSincard(const CovContext& ctxt)
    : AKernel(ECov::SINCARD, ctxt)
  {
  }

  KernelSincard::KernelSincard(const KernelSincard& r)
    : AKernel(r)
  {
  }

  KernelSincard& KernelSincard::operator=(const KernelSincard& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelSincard::~KernelSincard() {}

  double KernelSincard::getScadef() const
  {
    return (20.371);
  }

  double KernelSincard::_evaluateCov(double h) const
  {
    static double MIN_SIN = 1.e-5;
    double cov = 1.;
    if (h > MIN_SIN) cov = sin(h) / h;
    return (cov);
  }

  String KernelSincard::getFormula() const
  {
    return "C(h)=\\frac{sin(\\frac{h}{a})}{\\frac{h}{a}}";
  }

  double KernelSincard::simulateTurningBand(
    double t0,
    TurningBandOperate& operTB) const
  {
    return operTB.cosineOne(t0);
  }
} // namespace gstlrn
