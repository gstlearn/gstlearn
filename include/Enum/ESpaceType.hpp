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

#define ENUM_SPACETYPE ESpaceType, RN, \
                     RN, 1,  "Euclidean Space", \
                     SN, 2,  "Geometry on Sphere"

ENUM_DECLARE(ENUM_SPACETYPE)
