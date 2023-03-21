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

#define ENUM_LAW ELaw, CONSTANT,\
                 CONSTANT,    0, "Constant",\
                 UNIFORM,     1, "Uniform",\
                 GAUSSIAN,    2, "Gaussian",\
                 EXPONENTIAL, 3, "Exponential",\
                 GAMMA,       4, "Gamma",\
                 STABLE,      5, "Stable",\
                 BETA1,       6, "Beta-1",\
                 BETA2,       7, "Beta-2"

ENUM_DECLARE(ENUM_LAW)

