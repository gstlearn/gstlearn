/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_GAUSSINV EGaussInv, EMP,\
                 EMP,    0, "Empirical",\
                 HMT,    1, "Using Hermite Polynomials",\
                 NN,     2, "Nearest Neighbor"

ENUM_DECLARE(ENUM_GAUSSINV)
