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
#include "Covariances/CovReg1D.hpp"

#include "Covariances/CovContext.hpp"

CovReg1D::CovReg1D(const CovContext& ctxt)
: ACovFunc(ECov::REG1D, ctxt)
{
}

CovReg1D::CovReg1D(const CovReg1D &r)
: ACovFunc(r)
{
}

CovReg1D& CovReg1D::operator=(const CovReg1D &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovReg1D::~CovReg1D()
{
}

double CovReg1D::getScadef() const
{
  return (2.);
}

double CovReg1D::_evaluateCov(double h) const
{
  double cov = 0.;
  if (h < 1.)
  {
    cov = 1. - 3. * h * (1. - h / 2. * (1. + h / 6.));
  }
  else if (h < 2.)
  {
    cov = -2. + 3. * h * (1. - h / 2. * (1. - h / 6.));
  }
  return (cov);
}

