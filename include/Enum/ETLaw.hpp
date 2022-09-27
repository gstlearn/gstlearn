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

#define ENUM_TLAW ETLaw, CONSTANT,\
                 CONSTANT,    0, "Constant",\
                 UNIFORM,     1, "Uniform",\
                 GAUSSIAN,    2, "Gaussian",\
                 EXPONENTIAL, 3, "Exponential",\
                 GAMMA,       4, "Gamma",\
                 STABLE,      5, "Stable",\
                 BETA1,       6, "Beta-1",\
                 BETA2,       7, "Beta-2"

ENUM_DECLARE(ENUM_TLAW)

