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
#include "Covariances/CovLinear.hpp"

#include "Covariances/CovContext.hpp"
#include "Simulation/TurningBandOperate.hpp"

CovLinear::CovLinear(const CovContext& ctxt)
: ACovFunc(ECov::LINEAR, ctxt)
{
}

CovLinear::CovLinear(const CovLinear &r)
: ACovFunc(r)
{
}

CovLinear& CovLinear::operator=(const CovLinear &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovLinear::~CovLinear()
{
}

double CovLinear::_evaluateCov(double h) const
{
  double cov, r;

  r = getContext().getField();

  int ndim = getContext().getNDim();
  if (ndim == 1)
    cov = r - h;
  else if (ndim == 2)
    cov = r * GV_PI / 2 - h;
  else
    cov = r * 2 - h;

  return (cov);
}

double CovLinear::simulateTurningBand(double t0, TurningBandOperate &operTB) const
{
  return operTB.IRFProcessOne(t0);
}
