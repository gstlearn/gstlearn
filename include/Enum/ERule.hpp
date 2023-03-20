/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_RULE ERule, STD, \
                  STD,    0,  "Standard Lithotype Rule", \
                  SHIFT,  1,  "Shift Lithotype Rule", \
                  SHADOW, 2,  "Shadow Lithotype Rule"

ENUM_DECLARE(ENUM_RULE)
