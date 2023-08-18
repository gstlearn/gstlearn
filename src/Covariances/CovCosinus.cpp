/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovCosinus.hpp"

#include <math.h>
#include "Covariances/CovContext.hpp"

CovCosinus::CovCosinus(const CovContext& ctxt)
: ACovFunc(ECov::COSINUS, ctxt)
{
}

CovCosinus::CovCosinus(const CovCosinus &r)
: ACovFunc(r)
{
}

CovCosinus& CovCosinus::operator=(const CovCosinus &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovCosinus::~CovCosinus()
{
}

double CovCosinus::_evaluateCov(double h) const
{
  double cov = cos(2. * GV_PI * h);
  return (cov);
}

