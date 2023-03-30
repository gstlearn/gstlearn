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

class GSTLEARN_EXPORT IOptimCost {

public:
  IOptimCost() {};
  virtual ~IOptimCost() {};
  virtual void calculateGradient(const VectorDouble& indic,
                                 const VectorDouble& sval,
                                 double* normgrad) = 0;
};
