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
#include "Covariances/CovTriangle.hpp"

#include "Covariances/CovContext.hpp"

CovTriangle::CovTriangle(const CovContext& ctxt)
: ACovFunc(ECov::TRIANGLE, ctxt)
{
}

CovTriangle::CovTriangle(const CovTriangle &r)
: ACovFunc(r)
{
}

CovTriangle& CovTriangle::operator=(const CovTriangle &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovTriangle::~CovTriangle()
{
}

double CovTriangle::_evaluateCov(double h) const
{
  double cov = MAX(0, 1. - h);
  return (cov);
}

