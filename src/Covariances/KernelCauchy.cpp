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
#include "Covariances/KernelCauchy.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
KernelCauchy::KernelCauchy(const CovContext& ctxt)
  : AKernel(ECov::CAUCHY, ctxt)
{
  setParam(1);
}

KernelCauchy::KernelCauchy(const KernelCauchy& r)
  : AKernel(r)
{
}

KernelCauchy& KernelCauchy::operator=(const KernelCauchy& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelCauchy::~KernelCauchy()
{
}

double KernelCauchy::getScadef() const
{
  return sqrt(pow(20., 1. / getParam()) - 1.);
}

double KernelCauchy::_evaluateCov(double h) const
{
  double cov = 1. / pow(1. + (h * h), getParam());
  return (cov);
}

String KernelCauchy::getFormula() const
{
  return "C(h)=\\frac{1}{\\left( 1+ \\frac{h^2}{a_t^2} \\right)^\\alpha";
}
} // namespace gstlrn