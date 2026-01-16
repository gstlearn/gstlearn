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

#define ENUM_POST_UPSCALE EPostUpscale, MEAN, \
                 UNKNOWN, -1, "Unknown", \
                 NUM,      0, "Counter", \
                 MEAN,     1, "Average", \
                 MINI,     4, "Minimum", \
                 MAXI,     5, "Maximum"

ENUM_DECLARE(ENUM_POST_UPSCALE)
