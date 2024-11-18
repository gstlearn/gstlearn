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

#define ENUM_CONS_ELEM EConsElem, UNKNOWN, \
                       UNKNOWN,  0, "Unknown constraint", \
                       RANGE,    1, "Range", \
                       ANGLE,    2, "Anisotropy rotation angle (degree)", \
                       PARAM,    3, "Auxiliary parameter", \
                       SILL,     4, "Sill", \
                       SCALE,    5, "Scale", \
                       T_RANGE,  6, "Tapering range", \
                       VELOCITY, 7, "Velocity (advection)", \
                       SPHEROT,  8, "Rotation angle for Sphere", \
                       TENSOR,   9, "Anisotropy Matrix term"

ENUM_DECLARE(ENUM_CONS_ELEM)
