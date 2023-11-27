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
#include "Covariances/CovGamma.hpp"

#include "math.h"
#include "Covariances/CovContext.hpp"

CovGamma::CovGamma(const CovContext& ctxt)
: ACovFunc(ECov::GAMMA, ctxt)
{
  setParam(1);
}

CovGamma::CovGamma(const CovGamma &r)
: ACovFunc(r)
{
}

CovGamma& CovGamma::operator=(const CovGamma &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovGamma::~CovGamma()
{
}

double CovGamma::getScadef() const
{
  double param = getParam();
  if (param < 0.05) param = 1.; // This test is too avoid passing too small number to next line
  double scadef = pow(20.,1. / param) - 1.;
  return (scadef);
}

double CovGamma::_evaluateCov(double h) const
{
  double cov;
  cov = 1. / pow(1. + h, getParam());
  return (cov);
}

String CovGamma::getFormula() const
{
  return "C(h)=\\frac{1}{\\left( 1+ \\frac{h}{a_t} \\right)^\\alpha";
}
