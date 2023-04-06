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

#define ENUM_CONS_ELEM EConsElem, UNKNOWN, \
                       UNKNOWN,  0, "Unknown constraint", \
                       RANGE,    1, "Non-stationary range", \
                       ANGLE,    2, "Non-stationary anisotropy rotation angle (degree)", \
                       PARAM,    3, "Non-stationary auxiliary parameter", \
                       SILL,     4, "Non-stationary sill", \
                       SCALE,    5, "Non-stationary scale", \
                       T_RANGE,  6, "Non-stationary tapering range", \
                       VELOCITY, 7, "Non-stationary velocity (advection)", \
                       SPHEROT,  8, "Non-stationary rotation angle for Sphere", \
                       TENSOR,   9, "Non-stationary Anisotropy Matrix term"

ENUM_DECLARE(ENUM_CONS_ELEM)
