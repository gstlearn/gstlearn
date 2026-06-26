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
#define ENUM_POST_STAT                                                 \
  EPostStat, MEAN,                                                     \
  MEAN, 1, "Mean",                                                     \
  VAR, 2, "Variance",                                                  \
  VARP, 3, "Variance-P",                                               \
  STD, 4, "Std.",                                                      \
  STDP, 5, "Std-P",                                                    \
  MED, 6, "Median",                                                    \
  MINI, 7, "Minimum",                                                  \
  MAXI, 8, "Maximum"
// clang-format on

ENUM_DECLARE(ENUM_POST_STAT)
