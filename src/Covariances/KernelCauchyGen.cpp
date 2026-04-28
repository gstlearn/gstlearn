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
#include "Basic/Law.hpp"
#include "Basic/LawStable.hpp"
#include "Covariances/CovContext.hpp"
#include "Enum/ESimuType.hpp"

#include <cmath>

/*
C(h) = 1/(1 + abs(h)^alpha)^beta
with
alpha = _param[0] in (0,2]
beta  = _param[1] > 0
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

  KernelCauchyGen::~KernelCauchyGen() {}

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

  bool KernelCauchyGen::isValidForSimulation(const ESimuType& simuType) const
  {
    // spectral simulation is valid only for R^1
    return ((simuType == ESimuType::SPECTRAL) && (getContext().getNDim() == 1));
  }

  MatrixDense KernelCauchyGen::simulateSpectralOmega(Id nb) const
  {
    double alpha = getParam(0);
    double beta = getParam(1);
    MatrixDense omega(nb, 1);
    for (Id i = 0; i < nb; i++)
    {
      double value = pow(law_gamma(beta), 1 / alpha);
      if (alpha == 2.0)
      {
        value *= (sqrt(2.0) * law_gaussian());
      }
      else if (alpha <= 1.0)
      {
        value *= LawStable::law_cauchy();
        if (alpha < 1.0)
        {
          value *= LawStable::law_stable_unilateral(alpha);
        }
      }
      else if (alpha > 1.0)
      {
        value *= LawStable::law_stable_bilateral(alpha);
      }
      omega.setValue(i, 0, value);
    }
    return omega;
  }
} // namespace gstlrn
