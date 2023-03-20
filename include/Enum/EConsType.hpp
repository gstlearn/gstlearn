/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_CONS_TYPE EConsType, LOWER, \
                       LOWER,   -1, "Lower Bound", \
                       DEFAULT,  0, "Default parameter", \
                       UPPER,    1, "Upper Bound", \
                       EQUAL,    2, "Equality"

ENUM_DECLARE(ENUM_CONS_TYPE)
