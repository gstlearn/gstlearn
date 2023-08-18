/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovPenta.hpp"

#include "Covariances/CovContext.hpp"

CovPenta::CovPenta(const CovContext& ctxt)
: ACovFunc(ECov::PENTA, ctxt)
{
}

CovPenta::CovPenta(const CovPenta &r)
: ACovFunc(r)
{
}

CovPenta& CovPenta::operator=(const CovPenta &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovPenta::~CovPenta()
{
}

double CovPenta::_evaluateCov(double h) const
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

