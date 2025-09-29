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
#include "geoslib_define.h"

namespace gstlrn
{
GSTLEARN_EXPORT Id fftn(Id ndim,
                         const Id dims[],
                         double Re[],
                         double Im[],
                         Id iSign = 1,
                         double scaling = 1.);
}
