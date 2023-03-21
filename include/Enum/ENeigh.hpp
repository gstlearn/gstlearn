/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_NEIGH ENeigh, UNIQUE, \
                   UNIQUE, 0, "Unique Neighborhood", \
                   BENCH,  1, "Bench Neighborhood", \
                   MOVING, 2, "Moving Neighborhood", \
                   IMAGE,  3, "Image Neighborhood"

ENUM_DECLARE(ENUM_NEIGH)
