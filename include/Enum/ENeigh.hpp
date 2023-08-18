/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_NEIGH ENeigh, UNIQUE, \
                   UNKNOWN, -1, "Unknown Neighborhood", \
                   UNIQUE,   0, "Unique Neighborhood", \
                   BENCH,    1, "Bench Neighborhood", \
                   MOVING,   2, "Moving Neighborhood", \
                   CELL,     3, "Cell Neighborhood", \
                   IMAGE,    4, "Image Neighborhood"

ENUM_DECLARE(ENUM_NEIGH)
