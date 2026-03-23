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
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT LinearOpHelper
  {
  public:
    static double powerIteration(ALinearOp* op, Id niter = 30);
  };

  class LH: public LinearOpHelper
  {
  };

} // namespace gstlrn
