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
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
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

  KernelSpherical::~KernelSpherical() {}

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
    return "C(h)=1-\\frac{3}{2}\\left(\\frac{h}{a}\\right)+ "
           "\\frac{1}{2}\\left(\\frac{h}{a}\\right)^3";
  }

  double KernelSpherical::simulateTurningBand(
    double t0,
    TurningBandOperate& operTB) const
  {
    return operTB.shotNoiseAffineOne(t0);
  }

  double KernelSpherical::evaluateSpectrum(double freq) const
  {
    Id ndim = getContext().getNDim();

    // Handle zero frequency to avoid division by zero
    if (freq <= 0.0) return 1.0;

    double x = freq;
    // General order for isotropic Fourier transform in ndim dimensions
    double order = (static_cast<double>(ndim) - 2.0) / 2.0;

    // Compute cylindrical Bessel function J_order(x) (requires C++17 <cmath>)
    double bessel_val = cylBesselJ(order, x);

    // Normalized spectral density scaling according to dimension
    double val = pow(bessel_val / pow(x, order), 2.0);

    return val;
  }

  MatrixDense KernelSpherical::simulateSpectralOmega(Id nb) const
  {
    auto ndim = static_cast<Id>(getContext().getNDim());
    MatrixDense mat(nb, ndim);

    // Adjusted scaling factor to stretch the range to 12-13 while preserving origin slope
    double scaling_correction = 13.5;

    Id irow = 0;
    while (irow < nb)
    {
      // Draw raw frequency with a uniform proposal domain
      double raw_freq = law_uniform(0.0, 80.0);
      double target_density = evaluateSpectrum(raw_freq);

      // Algebraic proposal density matching the spherical spectrum tail
      double prop_density = 1.0 / (1.0 + raw_freq * raw_freq);
      double c_envelope = evaluateSpectrum(0.0) * 18.0;

      double u = law_uniform(0.0, 1.0);

      // Acceptance-rejection test
      if (u * (c_envelope * prop_density) <= target_density)
      {
        // Scale down the frequency using the adjusted correction factor
        double freq = raw_freq / scaling_correction;

        // Draw a uniform random direction on the unit sphere in ndim dimensions
        double norm_sq = 0.0;
        std::vector<double> dir(ndim);
        for (Id icol = 0; icol < ndim; icol++)
        {
          dir[icol] = law_gaussian();
          norm_sq += dir[icol] * dir[icol];
        }
        double norm = sqrt(norm_sq);

        // Assign scaled frequency vector components
        if (norm > 0)
        {
          for (Id icol = 0; icol < ndim; icol++)
          {
            double val = (dir[icol] / norm) * freq;
            mat.setValue(irow, icol, val);
          }
          irow++;
        }
      }
    }
    return mat;
  }
} // namespace gstlrn
