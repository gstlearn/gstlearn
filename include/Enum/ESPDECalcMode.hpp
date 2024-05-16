/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_SPDE_CALC_MODE ESPDECalcMode, KRIGING, \
                            KRIGING,      0,  "Kriging or Likelihood", \
                            KRIGVAR,      1,  "Kriging and St. Dev;", \
                            SIMUCOND,     2,  "Conditional simulations", \
                            SIMUNONCOND,  3,  "Non conditional simulations"

ENUM_DECLARE(ENUM_SPDE_CALC_MODE)
