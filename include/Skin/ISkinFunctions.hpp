/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

class GSTLEARN_EXPORT ISkinFunctions
{
public:
  ISkinFunctions() {};
  virtual ~ISkinFunctions() {};

  virtual int isAlreadyFilled(int /*ipos*/) const = 0;
  virtual int isToBeFilled(int /*ipos*/) const = 0;
  virtual double getWeight(int /*ipos*/, int /*idir*/) const { return 1; }
};
