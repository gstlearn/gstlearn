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
#include "LinearOp/ALinearOp.hpp"
#include "LinearOp/ASimulable.hpp"
#include "gstlearn_export.hpp"

class GSTLEARN_EXPORT ISimulable : public ALinearOp, public ASimulable
{
  public:
  ISimulable() { }
  virtual ~ISimulable() { }

};