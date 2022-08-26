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

#define ENUM_ROTATION ERotation, SXYZ, \
    SXYZ, 0, "sxyz", \
    SXYX, 1, "sxyx", \
    SXZY, 2, "sxzy", \
    SXZX, 3, "sxzx", \
    SYZX, 4, "syzx", \
    SYZY, 5, "syzy", \
    SYXZ, 6, "syxz", \
    SYXY, 7, "syxy", \
    SZXY, 8, "szxy", \
    SZXZ, 9, "szxz", \
    SZYX, 10, "szyx", \
    SZYZ, 11, "szyz", \
    RZYX, 12, "rzyx", \
    RXYX, 13, "rxyx", \
    RYZX, 14, "ryzx", \
    RXZX, 15, "rxzx", \
    RXZY, 16, "rxzy", \
    RYZY, 17, "ryzy", \
    RZXY, 18, "rzxy", \
    RYXY, 19, "ryxy", \
    RYXZ, 20, "ryxz", \
    RZXZ, 21, "rzxz", \
    RXYZ, 22, "rxyz", \
    RZYZ, 23, "rzyz"

ENUM_DECLARE(ENUM_ROTATION)
