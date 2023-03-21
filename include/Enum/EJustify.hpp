/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_JUSTIFY EJustify, LEFT, \
                     LEFT,  -1,  "Left Justification", \
                     CENTER, 0,  "Center Justification", \
                     RIGHT,  1,  "Right Justification"

ENUM_DECLARE(ENUM_JUSTIFY)
