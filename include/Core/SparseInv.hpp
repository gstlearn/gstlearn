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
GSTLEARN_EXPORT Id sparseinv(Id n,
                              Id* Lp,
                              Id* Li,
                              double* Lx,
                              double* d,
                              Id* Up,
                              Id* Uj,
                              double* Ux,
                              Id* Zp,
                              Id* Zi,
                              double* Zx,
                              double* z,
                              Id* Zdiagp,
                              Id* Lmunch);
}
