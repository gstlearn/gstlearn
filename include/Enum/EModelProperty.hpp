/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_MODEL_PROPERTY EModelProperty, NONE,\
                            NONE, 0,  "No specific property (LMC)", \
                            CONV, 1,  "Convolution mode", \
                            ANAM, 2,  "Anamorphosis mode", \
                            TAPE, 3,  "Tapering mode", \
                            GRAD, 4,  "Gradient mode"

ENUM_DECLARE(ENUM_MODEL_PROPERTY)
