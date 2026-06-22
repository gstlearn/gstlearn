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

namespace gstlrn
{
  GSTLEARN_EXPORT Id fftn(
    Id ndim,
    const Id dims[],
    double Re[],
    double Im[],
    Id iSign = 1,
    double scaling = 1.);
} // namespace gstlrn
