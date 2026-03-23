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
#include <Eigen/IterativeLinearSolvers>

namespace gstlrn
{

  class GSTLEARN_EXPORT APreconditioner
  {
  public:
    APreconditioner()
      : _status(Eigen::Success)
    {
    }

    virtual ~APreconditioner() = default;
#ifndef SWIG
    Eigen::ComputationInfo info() const { return _status; }
#endif

  private:
    Eigen::ComputationInfo _status;
  };

} // namespace gstlrn
