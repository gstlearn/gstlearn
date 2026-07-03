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
#define ENUM_VARIABLENUMERIC                                            \
  EVariableNumeric, NONE,                                               \
  NONE,       0, "No specific Format",                                  \
  AUTO,       1, "Automatic Format",                                    \
  FIXED,      2, "Fixed Format",                                        \
  SCIENTIFIC, 3, "Scientific Format"
// clang-format on

ENUM_DECLARE(ENUM_VARIABLENUMERIC)
