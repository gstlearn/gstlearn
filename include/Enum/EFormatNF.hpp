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
#define ENUM_FORMATNF                                                  \
  EFormatNF, DEFAULT,                                                  \
  DEFAULT, 0, "Using the Defaulted Format",                            \
  ASCII, 1, "Using ASCII Format",                                      \
  H5, 2, "Using the HDF5 Format"
// clang-format on

ENUM_DECLARE(ENUM_FORMATNF)
