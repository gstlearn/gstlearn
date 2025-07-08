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
#include "Covariances/CovGaussian.hpp"
#include "Basic/Law.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Simulation/TurningBandOperate.hpp"

#include <cmath>

namespace gstlrn
{
CovGaussian::CovGaussian(const CovContext& ctxt)
  : ACovFunc(ECov::GAUSSIAN, ctxt)
{
}

CovGaussian::CovGaussian(const CovGaussian& r)
  : ACovFunc(r)
{
}

CovGaussian& CovGaussian::operator=(const CovGaussian& r)
{
  if (this != &r)
  {
    ACovFunc::operator=(r);
  }
  return *this;
}

CovGaussian::~CovGaussian()
{
}

double CovGaussian::getScadef() const
{
  return (1.730818);
}

double CovGaussian::_evaluateCov(double h) const
{
  double r2 = h * h;
  if (r2 > MAX_EXP) return 0.;
  double cov = exp(-r2);
  return (cov);
}

double CovGaussian::_evaluateCovDerivative(int degree, double h) const
{
  double h2 = h * h;
  if (h2 > MAX_EXP) return 0.;

  double cov = 0.;
  switch (degree)
  {
    case 1: // First order derivative
      cov = -2. * h * exp(-h2); // Derivative of e^(-h^2) is -2h * e^(-h^2)
      break;

    case 2: // Second order derivative
      cov = (4. * h2 - 2.) * exp(-h2); // 
      break;

    case 3: // Third order derivative
      cov = 4. * exp(-h2) * h * (3 - 2. * h2); //[-2h * (4h^2-2) +8h]e^(-h^2) = (12h -8h^3)e(-h2)  
      break;

    case 4: // Fourth order derivative
      double h4 = h2 * h2;
      cov       = 4. * exp(-h2) * (3. - 12. * h2 + 4 * h4);  //[(12-24h2) -2h(12h-8h^3)]e(-h^2) = (12-48h2 + 16h4)
      break;
  }
  return (cov);
}

String CovGaussian::getFormula() const
{
  return "C(h)=1-e^{-h^2}";
}

double CovGaussian::simulateTurningBand(double t0, TurningBandOperate& operTB) const
{
  return operTB.cosineOne(t0);
}

MatrixDense CovGaussian::simulateSpectralOmega(int nb) const
{
  int ndim = getContext().getNDim();
  MatrixDense mat(nb, ndim);
  for (int irow = 0; irow < nb; irow++)
    for (int icol = 0; icol < ndim; icol++)
      mat.setValue(irow, icol, law_gaussian() * sqrt(2.));
  return mat;
}
} // namespace gstlrn