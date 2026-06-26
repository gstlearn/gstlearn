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
#define ENUM_DIRGEN                                                    \
  EDirGen, VDC,                                                        \
  VDC, 0, "Van der Corput",                                            \
  RND, 1, "Random Directions"
// clang-format on

ENUM_DECLARE(ENUM_DIRGEN)
