/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_KRIG_OPT EKrigOpt, POINT, \
                      POINT,    0,  "Punctual estimation", \
                      BLOCK,    1,  "Block average estimation", \
                      DRIFT,    2,  "Large scale Drift estimation", \
                      DGM,      3,  "Discrete Gaussian Model"

ENUM_DECLARE(ENUM_KRIG_OPT)
