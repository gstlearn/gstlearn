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
#include "Covariances/KernelNugget.hpp"
#include "Covariances/CovContext.hpp"

namespace gstlrn
{
KernelNugget::KernelNugget(const CovContext& ctxt)
  : AKernel(ECov::NUGGET, ctxt)
{
}

KernelNugget::KernelNugget(const KernelNugget& r)
  : AKernel(r)
{
}

KernelNugget& KernelNugget::operator=(const KernelNugget& r)
{
  if (this != &r)
  {
    AKernel::operator=(r);
  }
  return *this;
}

KernelNugget::~KernelNugget()
{
}

double KernelNugget::_evaluateCov(double h) const
{
  double cov = 0.;
  if (ABS(h) < 1.e-10) cov = 1.;
  return (cov);
}

String KernelNugget::getFormula() const
{
  return "C(h)=\\delta(h)";
}
}