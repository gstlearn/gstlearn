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
#include "LinearOp/ASimulable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "geoslib_define.h"

VectorDouble ASimulable::evalSimulate(const VectorDouble& whitenoise) const
{
  VectorDouble res;
  evalSimulate(whitenoise, res);
  return res;
}

VectorDouble ASimulable::simulate() const
{
  VectorDouble whitenoise(getSize());
  VH::simulateGaussianInPlace(whitenoise);
  VectorDouble res(getSize());
  evalSimulate(whitenoise, res);
  return res;
}
int ASimulable::evalSimulate(const VectorDouble& whitenoise,
                             VectorDouble& outv) const
{
  outv.resize(whitenoise.size());
  constvect ws(whitenoise);
  vect outs(outv);
  return evalSimulate(ws, outs);
}

int ASimulable::evalSimulate(const constvect whitenoise, vect result) const
{
  std::fill(result.begin(),result.end(),0.);
  return _addSimulateToDest(whitenoise, result);
}

int ASimulable::addSimulateToDest(const constvect whitenoise, vect outv) const
{
  return _addSimulateToDest(whitenoise, outv);
}

double ASimulable::computeLogDet(int nMC) const
{
  DECLARE_UNUSED(nMC);
  messerr("computeLogDet not implemented in ASimulable"); 
  return TEST;
}
