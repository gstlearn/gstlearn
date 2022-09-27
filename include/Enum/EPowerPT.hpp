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

#define ENUM_POWER_PT EPowerPT, UNDEFINED, \
                      UNDEFINED, -1,  "Power is undefined", \
                      ONE,        0,  "Power is 1", \
                      MINUSONE,   1,  "Power is -1", \
                      MINUSHALF,  2,  "Power is -0.5", \
                      HALF,       3,  "Power is 0.5", \
                      LOG,        4,  "Logarithm"

ENUM_DECLARE(ENUM_POWER_PT)
