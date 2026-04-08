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

#define ENUM_CSV                                                               \
  ECSV, ENGLISH, FRENCH, 0, "French CSV", ENGLISH, 1, "English CSV",           \
    TABULATED, 2, "Tabulated CSV"

ENUM_DECLARE(ENUM_CSV)
