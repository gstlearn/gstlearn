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
#define ENUM_SELECTIVITY                                               \
  ESelectivity, UNDEFINED,                                             \
  UNDEFINED, -1, "Unknown Option",                                     \
  Z, 0, "Grade",                                                       \
  T, 1, "Tonnage",                                                     \
  Q, 2, "Metal quantity",                                              \
  B, 3, "Conventional Benefit",                                        \
  M, 4, "Recovered mean",                                              \
  PROP, 5, "Probability to exceed Cutoff",                             \
  QUANT, 6, "Quantile"
// clang-format on

ENUM_DECLARE(ENUM_SELECTIVITY)
