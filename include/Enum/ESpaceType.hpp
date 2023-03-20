/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_SPACETYPE ESpaceType, RN, \
                     RN, 1,  "Euclidean Space", \
                     SN, 2,  "Geometry on Sphere"

ENUM_DECLARE(ENUM_SPACETYPE)
