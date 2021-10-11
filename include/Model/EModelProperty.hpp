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

#define ENUM_MODEL_PROPERTY EModelProperty, NONE,\
                            NONE, 0,  "No specific property", \
                            CONV, 1,  "Convolution mode", \
                            ANAM, 2,  "Anamorphosis mode", \
                            TAPE, 3,  "Tapering mode", \
                            GRAD, 4,  "Gradient mode"

ENUM_DECLARE(ENUM_MODEL_PROPERTY)
