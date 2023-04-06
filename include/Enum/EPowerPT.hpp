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

#define ENUM_POWER_PT EPowerPT, UNDEFINED, \
                      UNDEFINED, -1,  "Power is undefined", \
                      ONE,        0,  "Power is 1", \
                      MINUSONE,   1,  "Power is -1", \
                      MINUSHALF,  2,  "Power is -0.5", \
                      HALF,       3,  "Power is 0.5", \
                      LOG,        4,  "Logarithm"

ENUM_DECLARE(ENUM_POWER_PT)
