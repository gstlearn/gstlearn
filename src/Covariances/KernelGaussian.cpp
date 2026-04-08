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
#include "Covariances/KernelGaussian.hpp"
#include "Basic/Law.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Simulation/TurningBandOperate.hpp"

#include <cmath>
#include <cstddef>

namespace gstlrn
{
  KernelGaussian::KernelGaussian(const CovContext& ctxt)
    : AKernel(ECov::GAUSSIAN, ctxt)
  {
  }

  KernelGaussian::KernelGaussian(const KernelGaussian& r)
    : AKernel(r)
  {
  }

  KernelGaussian& KernelGaussian::operator=(const KernelGaussian& r)
  {
    if (this != &r)
    {
      AKernel::operator=(r);
    }
    return *this;
  }

  KernelGaussian::~KernelGaussian() {}

  double KernelGaussian::getScadef() const
  {
    return (1.730818);
  }

  double KernelGaussian::_evaluateCov(double h) const
  {
    double r2 = h * h;
    if (r2 > MAX_EXP) return 0.;
    double cov = exp(-r2);
    return (cov);
  }

  double KernelGaussian::_evaluateCovFirstDerivativeOverH(double h) const
  {
    double r2 = h * h;
    if (r2 > MAX_EXP) return 0.;
    double cov = -2. * exp(-r2);
    return cov;
  }

  double KernelGaussian::_evaluateCovDerivative(Id degree, double h) const
  {
    double h2 = h * h;
    if (h2 > MAX_EXP) return 0.;

    double cov = 0.;
    double expr2 = exp(-h2);
    switch (degree)
    {
      case 1: // First order derivative
        cov = -2. * h * expr2; // Derivative of e^(-h^2) is -2h * e^(-h^2)
        break;

      case 2: // Second order derivative
        cov =
          (4. * h2 - 2.) * expr2; // Derivative of e^(-h^2) is -2h * e^(-h^2)
        break;

      case 3: // Third order derivative
        cov = 4. * expr2 * h
            * (3 - 2. * h2); //[-2h * (4h^2-2) +8h]e^(-h^2) = (12h -8h^3)e(-h2)
        break;

      case 4: // Fourth order derivative
        double h4 = h2 * h2;
        cov = 4. * expr2
            * (3. - 12. * h2
               + 4 * h4); //[(12-24h2) -2h(12h-8h^3)]e(-h^2) = (12-48h2 + 16h4)
        break;
    }
    return (cov);
  }

  String KernelGaussian::getFormula() const
  {
    return "C(h)=e^{-h^2}";
  }

  double KernelGaussian::simulateTurningBand(
    double t0,
    TurningBandOperate& operTB) const
  {
    return operTB.cosineOne(t0);
  }

  double KernelGaussian::evaluateSpectrum(double freq) const
  {
    size_t ndim = getContext().getNDim();
    double val = exp(-freq * freq * 0.25) / pow(2 * sqrt(GV_PI), ndim);
    return val;
  }

  MatrixDense KernelGaussian::simulateSpectralOmega(Id nb) const
  {
    size_t ndim = getContext().getNDim();
    MatrixDense mat(nb, ndim);
    double sqrt2 = sqrt(2.0);
    for (size_t icol = 0; icol < ndim; icol++)
    {
      auto view = mat.getViewOnColumnModify(icol);
      for (size_t irow = 0; irow < static_cast<size_t>(nb); irow++)
        view[irow] = sqrt2 * law_gaussian();
    }
    return mat;
  }
} // namespace gstlrn
