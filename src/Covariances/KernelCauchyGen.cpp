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

namespace gstlrn
{
KernelCauchyGen::KernelCauchyGen(const CovContext& ctxt)
  : AKernel(ECov::CAUCHY, ctxt)
{
  setParam(1);
  setAlpha(2.0);
}

KernelCauchyGen::KernelCauchyGen(const KernelCauchyGen& r)
  : AKernel(r)
  , _alpha(r._alpha)
{
}

KernelCauchyGen& KernelCauchyGen::operator=(const KernelCauchyGen& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
    _alpha = r._alpha;
  }
  return *this;
}

KernelCauchyGen::~KernelCauchyGen()
{
}

bool KernelCauchyGen::setAlpha(double alpha)
{
    if((alpha > 0)&&(alpha <= 2.0)){
        _alpha = alpha;
        return true;
    }
  _alpha = 0.0;
  return false;
}

double KernelCauchyGen::getScadef() const
{
  return pow(pow(20., 1. / getParam()) - 1., 1. / getAlpha());
}

double KernelCauchyGen::_evaluateCov(double h) const
{
  double cov = 1. / pow(1. + pow(abs(h), getAlpha()), getParam());
  return (cov);
}

String KernelCauchyGen::getFormula() const
{
  return "C(h)=\\frac{1}{\\left( 1+ \\frac{h^\\alpha}{a_t^\\alpha} \\right)^\\beta";
}
} // namespace gstlrn