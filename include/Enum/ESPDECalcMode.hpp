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

#define ENUM_SPDE_CALC_MODE ESPDECalcMode, KRIGING, \
                            KRIGING,      0,  "Kriging", \
                            SIMUCOND,     1,  "Conditional simulations", \
                            SIMUNONCOND,  2,  "Non conditional simulations", \
                            LIKELIHOOD,   3,  "Likelihood computations"

ENUM_DECLARE(ENUM_SPDE_CALC_MODE)
