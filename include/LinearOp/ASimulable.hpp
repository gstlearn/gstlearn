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
#pragma once
#include "gstlearn_export.hpp"
#include "Basic/VectorNumT.hpp"
#include "LinearOp/ALinearOp.hpp"

namespace gstlrn{
class GSTLEARN_EXPORT ASimulable : public virtual ALinearOp
{
public:
  ASimulable();
  virtual ~ASimulable() {}

  Id evalSimulate(const VectorDouble& whitenoise, VectorDouble& outv) const;
  VectorDouble evalSimulate(const VectorDouble& whitenoise) const;
  VectorDouble simulate() const;
  virtual double computeLogDet(Id nMC = 1) const; 

#ifndef SWIG
public:
  Id evalSimulate(const constvect whitenoise, vect result) const;
  Id addSimulateToDest(const constvect whitenoise, vect outv) const;
protected:
  virtual Id _addSimulateToDest(const constvect whitenoise, vect outv) const = 0;
#endif
};
}