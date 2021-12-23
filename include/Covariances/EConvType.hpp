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

#define ENUM_CONVTYPE EConvType, UNKNOWN,\
                 UNKNOWN,    -1, "Unknown",\
                 UNIFORM,     0, "Uniform",\
                 EXPONENTIAL, 1, "Exponential",\
                 GAUSSIAN,    2, "Gaussian",\
                 SINCARD,     3, "Cardinal Sine"

ENUM_DECLARE(ENUM_CONVTYPE)
