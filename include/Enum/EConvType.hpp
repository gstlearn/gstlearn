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
#define ENUM_CONVTYPE                                                  \
  EConvType, UNIFORM,                                                  \
  UNIFORM, 0, "Uniform",                                               \
  EXPONENTIAL, 1, "Exponential",                                       \
  GAUSSIAN, 2, "Gaussian",                                             \
  SINCARD, 3, "Cardinal Sine"
// clang-format on

ENUM_DECLARE(ENUM_CONVTYPE)
