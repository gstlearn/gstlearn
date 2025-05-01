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
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Basic/Law.hpp"

#include "math.h"

CovGaussian::CovGaussian(const CovContext& ctxt)
: ACovFunc(ECov::GAUSSIAN, ctxt)
{
}

CovGaussian::CovGaussian(const CovGaussian &r)
: ACovFunc(r)
{
}

CovGaussian& CovGaussian::operator=(const CovGaussian &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
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
  double r2 = h * h;
  if (r2 > MAX_EXP) return 0.;

  double cov = 0.;
  switch (degree)
  {
    case 1:   // First order derivative
      cov = 2. * exp(-r2);
      break;

    case 2: // Second order derivative
      cov = (4. * r2 - 2.) * exp(-r2);
      break;

    case 3: // Third order derivative
      cov = 4. * exp(-r2) * h * (3 - 2. * r2);
      break;

    case 4: // Fourth order derivative
      double r4 = r2 * r2;
      cov = 8. * exp(-r2) * (6. - 15. * r2 + r4);
      break;
  }
  return (cov);
}

String CovGaussian::getFormula() const
{
  return "C(h)=1-\\frac{7h^2}{a^2}+\\frac{35h^3}{4a^3}-\\frac{7h^5}{2a^5}-\\frac{3h^7}{4a^7}";
}

double CovGaussian::simulateTurningBand(double t0, TurningBandOperate &operTB) const
{
  return operTB.cosineOne(t0);
}

MatrixDense CovGaussian::simulateSpectralOmega(int nb) const
{
  int ndim = getContext().getNDim();
  MatrixDense mat(nb, ndim);
  for (int irow = 0; irow < nb; irow++)
    for (int icol = 0; icol < ndim; icol++)
      mat.setValue(irow, icol, law_gaussian());
  return mat;
}
