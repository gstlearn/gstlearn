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
#include "Covariances/CovWendland1.hpp"
#include "Covariances/CovContext.hpp"

namespace gstlrn
{
CovWendland1::CovWendland1(const CovContext& ctxt)
  : ACovFunc(ECov::WENDLAND1, ctxt)
{
}

CovWendland1::CovWendland1(const CovWendland1& r)
  : ACovFunc(r)
{
}

CovWendland1& CovWendland1::operator=(const CovWendland1& r)
{
  if (this != &r)
  {
    ACovFunc::operator=(r);
  }
  return *this;
}

CovWendland1::~CovWendland1()
{
}

double CovWendland1::_evaluateCov(double h) const
{
  // From "Computed Supported Correlation Functions" by T. Gneiting with n=3
  double cov = 0.;
  if (h < 1)
    cov = (1. + 4. * h) * pow(1. - h, 4.);
  return (cov);
}

}