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
#include "Covariances/KernelCauchyGen.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

/*
C(h) = 1/(1 + abs(h)^alpha)^beta
with
alpha = _param[0]
beta  = _param[1]
*/
namespace gstlrn
{
KernelCauchyGen::KernelCauchyGen(const CovContext& ctxt)
  : AKernel(ECov::CAUCHY_GEN, ctxt)
{
  setParam(2, 0);
  setParam(1, 1);
}

KernelCauchyGen::KernelCauchyGen(const KernelCauchyGen& r)
  : AKernel(r)
{
}

KernelCauchyGen& KernelCauchyGen::operator=(const KernelCauchyGen& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelCauchyGen::~KernelCauchyGen()
{
}

double KernelCauchyGen::getScadef() const
{
  return pow(pow(20., 1. / getParam(1)) - 1., 1. / getParam(0));
}

double KernelCauchyGen::_evaluateCov(double h) const
{
  double cov = 1. / pow(1. + pow(abs(h), getParam(0)), getParam(1));
  return (cov);
}

String KernelCauchyGen::getFormula() const
{
  return "C(h)=\\frac{1}{\\left( 1+ h^\\alpha \\right)^\\beta";
}
} // namespace gstlrn