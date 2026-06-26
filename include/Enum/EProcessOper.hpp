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
#define ENUM_PROCESS_OPER                                              \
  EProcessOper, COPY,                                                  \
  UNDEFINED, -1, "Undefined Process Operation",                        \
  COPY, 0, "Copy Process",                                             \
  MARGINAL, 1, "Marginal Process",                                     \
  CONDITIONAL, 2, "Conditional Process"
// clang-format on

ENUM_DECLARE(ENUM_PROCESS_OPER)
