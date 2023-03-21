/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
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
