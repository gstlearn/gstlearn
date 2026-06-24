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
#define ENUM_SIMUTYPE                                                  \
  ESimuType, UNDEFINED,                                                \
  UNDEFINED, 0, "Not defined yet",                                     \
  TB, 1, "Turning Bands method",                                       \
  SPECTRAL, 2, "Spectral method"
// clang-format on

ENUM_DECLARE(ENUM_SIMUTYPE)
