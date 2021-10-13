/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
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
                       SPHEROT,  8, "Non-stationary rotation angle for Sphere"

ENUM_DECLARE(ENUM_CONS_ELEM)
