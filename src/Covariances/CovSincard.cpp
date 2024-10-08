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
#include "Covariances/CovSincard.hpp"
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"
#include "Basic/Law.hpp"

#include "math.h"

CovSincard::CovSincard(const CovContext& ctxt)
: ACovFunc(ECov::SINCARD, ctxt)
{
}

CovSincard::CovSincard(const CovSincard &r)
: ACovFunc(r)
{
}

CovSincard& CovSincard::operator=(const CovSincard &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovSincard::~CovSincard()
{
}

double CovSincard::getScadef() const
{
  return (20.371);
}

double CovSincard::_evaluateCov(double h) const
{
  static double MIN_SIN = 1.e-5;
  double cov = 1.;
  if (h > MIN_SIN) cov = sin(h) / h;
  return (cov);
}

String CovSincard::getFormula() const
{
  return "C(h)=\\frac{sin(\\frac{h}{a})}{\\frac{h}{a}}";
}

double CovSincard::simulateTurningBand(double t0, TurningBandOperate &operTB) const
{
  return operTB.cosineOne(t0);
}
