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
#define ENUM_MODEL_PROPERTY                                            \
  EModelProperty, CONV,                                                \
  NONE, 0, "No specific property (LMC)",                               \
  CONV, 1, "Convolution mode",                                         \
  ANAM, 2, "Anamorphosis mode",                                        \
  TAPE, 3, "Tapering mode"
// clang-format on

ENUM_DECLARE(ENUM_MODEL_PROPERTY)
