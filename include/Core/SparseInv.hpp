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

GSTLEARN_EXPORT int sparseinv(int n,
                              int* Lp,
                              int* Li,
                              double* Lx,
                              double* d,
                              int* Up,
                              int* Uj,
                              double* Ux,
                              int* Zp,
                              int* Zi,
                              double* Zx,
                              double* z,
                              int* Zdiagp,
                              int* Lmunch);
