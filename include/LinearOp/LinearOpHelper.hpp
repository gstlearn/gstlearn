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
#include "LinearOp/ALinearOp.hpp"

namespace gstlrn
{
class GSTLEARN_EXPORT LinearOpHelper
{
  public:
    static double powerIteration(const ALinearOp* op, Id niter = 30);
    
};

class LH : public LinearOpHelper
{

};

}