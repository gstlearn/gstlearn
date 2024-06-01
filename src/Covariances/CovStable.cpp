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
#include "Covariances/CovStable.hpp"

#include <math.h>
#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"
#include "Basic/Law.hpp"

CovStable::CovStable(const CovContext& ctxt)
: ACovFunc(ECov::STABLE, ctxt)
{
  setParam(1);
}

CovStable::CovStable(const CovStable &r)
: ACovFunc(r)
{
}

CovStable& CovStable::operator=(const CovStable &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovStable::~CovStable()
{
}

double CovStable::getScadef() const
{
  return pow(3., 1. / getParam());
}

double CovStable::_evaluateCov(double h) const
{
  double cov = 1.;
  if (h > 0) cov = exp(-pow(h, getParam()));

  return (cov);
}

double CovStable::simulateTurningBand(double t0,
                                      const VectorDouble &t,
                                      TurningBandOperate &operTB) const
{
  if (getParam() > 1)
    return operTB.cosineOne(t0, t);
  else
    return operTB.spectralOne(t0, t);
}
