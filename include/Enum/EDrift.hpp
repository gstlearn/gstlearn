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

#define ENUM_DRIFT EDrift, UNKNOWN, \
                   UNKNOWN, -1, "Unknown Drift", \
                   UC,       0, "Universality Condition", \
                   X,        1, "Drift along X", \
                   Y,        2, "Drift along Y", \
                   Z,        3, "Drift along Z", \
                   X2,       4, "Drift along X", \
                   Y2,       5, "Drift along Y²", \
                   XY,       6, "Drift along XY", \
                   Z2,       7, "Drift along Z²", \
                   XZ,       8, "Drift along XZ", \
                   YZ,       9, "Drift along YZ", \
                   X3,      10, "Drift along X³", \
                   X2Y,     11, "Drift along X²Y", \
                   XY2,     12, "Drift along XY²", \
                   Y3,      13, "Drift along Y³", \
                   Z3,      14, "Drift along Z³", \
                   F,       15, "External Drift"

ENUM_DECLARE(ENUM_DRIFT)
