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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "LinearOp/ASimulable.hpp"

namespace gstlrn
{
  class GSTLEARN_EXPORT IPrecisionOp: virtual public ASimulable
  {
  public:
    IPrecisionOp(); // default but see below
    IPrecisionOp(const IPrecisionOp&) = default;

    IPrecisionOp(IPrecisionOp&& m) noexcept
    {
      if (this != &m) ASimulable::operator=(std::move(m));
    }

    IPrecisionOp& operator=(const IPrecisionOp&) = default;

    IPrecisionOp& operator=(IPrecisionOp&& m) noexcept
    {
      if (this != &m) ASimulable::operator=(std::move(m));
      return *this;
    }

    virtual ~IPrecisionOp() = default;

    VectorDouble getRangeEigenVal(Id ndiscr = 100) const
    {
      std::pair<double, double> ranges = rangeEigenVal(ndiscr);
      VectorDouble result(2);
      result[0] = ranges.first;
      result[1] = ranges.second;
      return result;
    }

#ifndef SWIG
    virtual std::pair<double, double> rangeEigenVal(Id ndiscr = 100) const = 0;
#endif
  };

  // For some reason VS2017 crashes at runtime when a derived class calls this constructor
  // if = default is in the class declaration, but works if it's inline here (which makes
  // the constructor "user-provided" and may cause the difference?). This is a bug,
  // most likely due to a mix of dllexport, default keyword and virtual (base) classes?
  inline IPrecisionOp::IPrecisionOp() = default;

} // namespace gstlrn
