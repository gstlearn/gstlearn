/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_KRIG_OPT EKrigOpt, PONCTUAL, \
                      PONCTUAL, 0,  "Punctual estimation", \
                      BLOCK,    1,  "Block average estimation", \
                      DRIFT,    2,  "Large scale Drift estimation", \
                      DGM,      3,  "Discrete Gaussian Model"

ENUM_DECLARE(ENUM_KRIG_OPT)
