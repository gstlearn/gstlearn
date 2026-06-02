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

namespace gstlrn
{
  class GSTLEARN_EXPORT TracePGS
  {
  public:
    TracePGS(
      Id ngrf = 0,
      Id npar = 0,
      bool flag_rho = false,
      bool flag_correl = false);
    TracePGS(const TracePGS& r) = default;
    TracePGS& operator=(const TracePGS& r) = default;
    virtual ~TracePGS() = default;

    void define(Id ngrf, Id npar, Id flag_rho, bool flag_correl);
    double extractTrace();
    void addRow();
    void update(
      double value0,
      double value1,
      Id origin,
      Id number,
      const double* values);

  private:
    bool _flagTrace;
    Id _nrow;
    Id _ncol;
    VectorDouble _trace;
  };
} // namespace gstlrn
