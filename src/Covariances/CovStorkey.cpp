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
#include "Covariances/CovStorkey.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
CovStorkey::CovStorkey(const CovContext& ctxt)
  : ACovFunc(ECov::STORKEY, ctxt)
{
}

CovStorkey::CovStorkey(const CovStorkey& r)
  : ACovFunc(r)
{
}

CovStorkey& CovStorkey::operator=(const CovStorkey& r)
{
  if (this != &r)
  {
    ACovFunc::operator=(r);
  }
  return *this;
}

CovStorkey::~CovStorkey()
{
}

double CovStorkey::_evaluateCov(double h) const
{
  double cov = 0.;
  double pi2 = 2. * GV_PI;
  if (h < 1)
    cov = (2. * (1. - h) * (1. + cos(pi2 * h) / 2.) + 3 / pi2 * sin(pi2 * h)) / 3.;
  return (cov);
}

}