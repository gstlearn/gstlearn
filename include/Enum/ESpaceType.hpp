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
#define ENUM_SPACETYPE                                                 \
  ESpaceType, RN,                                                      \
  COMPOSITE, 0, "Composite Space",                                     \
  RN, 1, "Euclidean Space",                                            \
  SN, 2, "Geometry on Sphere"
// clang-format on

ENUM_DECLARE(ENUM_SPACETYPE)
