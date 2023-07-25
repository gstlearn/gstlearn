/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_POST_UPSCALE EPostUpscale, UNKNOWN, \
                 UNKNOWN, -1, "Unknown Option", \
                 NUM,      0, "Counter", \
                 MEAN,     1, "Average", \
                 MINI,     4, "Minimum", \
                 MAXI,     5, "Maximum"

ENUM_DECLARE(ENUM_POST_UPSCALE)
