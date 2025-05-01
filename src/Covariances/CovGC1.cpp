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
#include "Covariances/CovGC1.hpp"

#include "Simulation/TurningBandOperate.hpp"
#include "Covariances/CovContext.hpp"

CovGC1::CovGC1(const CovContext& ctxt)
: ACovFunc(ECov::ORDER1_GC, ctxt)
{
}

CovGC1::CovGC1(const CovGC1 &r)
: ACovFunc(r)
{
}

CovGC1& CovGC1::operator=(const CovGC1 &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovGC1::~CovGC1()
{
}

double CovGC1::_evaluateCov(double h) const
{
  double cov;
  double r = getContext().getField();
  int ndim = getContext().getNDim();

  if (ndim == 1)
    cov = r - h;
  else if (ndim == 2)
    cov = r * GV_PI / 2 - h;
  else
    cov = r * 2 - h;

  return (cov);
}

double CovGC1::simulateTurningBand(double t0, TurningBandOperate &operTB) const
{
  return operTB.IRFProcessOne(t0);
}
