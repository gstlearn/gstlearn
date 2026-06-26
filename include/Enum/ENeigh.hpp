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

#include "Enum/AEnum.hpp"

// clang-format off
#define ENUM_NEIGH                                                     \
  ENeigh, UNIQUE,                                                      \
  UNIQUE, 0, "Unique Neighborhood",                                    \
  BENCH, 1, "Bench Neighborhood",                                      \
  MOVING, 2, "Moving Neighborhood",                                    \
  CELL, 3, "Cell Neighborhood",                                        \
  IMAGE, 4, "Image Neighborhood"
// clang-format on

ENUM_DECLARE(ENUM_NEIGH)
