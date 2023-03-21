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

#define ENUM_SPACETYPE ESpaceType, RN, \
                     RN, 1,  "Euclidean Space", \
                     SN, 2,  "Geometry on Sphere"

ENUM_DECLARE(ENUM_SPACETYPE)
