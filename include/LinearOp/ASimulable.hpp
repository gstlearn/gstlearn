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

class GSTLEARN_EXPORT ASimulable : public ALinearOp
{
public:
  ASimulable() {}
  virtual ~ASimulable() {}

  int evalSimulate(const VectorDouble& whitenoise, VectorDouble& outv) const;
  VectorDouble evalSimulate(const VectorDouble& whitenoise) const;
  VectorDouble simulate() const;
#ifndef SWIG

public:
  int evalSimulate(const constvect whitenoise, vect result) const;
  int addSimulateToDest(const constvect whitenoise, vect outv) const;

protected:
  virtual int _addSimulateToDest(const constvect whitenoise, vect outv) const = 0;
#endif
};
