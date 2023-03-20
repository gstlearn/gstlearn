/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
