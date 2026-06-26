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
#define ENUM_CONS_TYPE                                                 \
  EConsType, LOWER,                                                    \
  LOWER, -1, "Lower Bound",                                            \
  DEFAULT, 0, "Default parameter",                                     \
  UPPER, 1, "Upper Bound",                                             \
  EQUAL, 2, "Equality"
// clang-format on

ENUM_DECLARE(ENUM_CONS_TYPE)
