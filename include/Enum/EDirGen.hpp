/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_DIRGEN EDirGen, VDC,\
                 VDC,    0, "Van der Corput",\
                 RND,    1, "Random Directions"

ENUM_DECLARE(ENUM_DIRGEN)
