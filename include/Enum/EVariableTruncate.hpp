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
#define ENUM_VARIABLETRUNCATE                                            \
  EVariableTruncate, RIGHT,                                              \
  LEFT,  0, "Align to the Left",                                         \
  RIGHT, 1, "Align to the Right"
// clang-format on

ENUM_DECLARE(ENUM_VARIABLETRUNCATE)
