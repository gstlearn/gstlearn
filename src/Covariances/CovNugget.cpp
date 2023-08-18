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
#include "Covariances/CovNugget.hpp"

#include "Covariances/CovContext.hpp"

CovNugget::CovNugget(const CovContext& ctxt)
: ACovFunc(ECov::NUGGET, ctxt)
{
}

CovNugget::CovNugget(const CovNugget &r)
: ACovFunc(r)
{
}

CovNugget& CovNugget::operator=(const CovNugget &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovNugget::~CovNugget()
{
}

double CovNugget::_evaluateCov(double h) const
{
  double cov = 0.;
  if (ABS(h) < 1.e-10) cov = 1.;
  return (cov);
}

String CovNugget::getFormula() const
{
  return "C(h)=\\delta(h)";
}
