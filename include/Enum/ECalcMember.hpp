/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_CALC_MEMBER ECalcMember, LHS, \
                             LHS, 0, "Left-hand Side of the Kriging System", \
                             RHS, 1, "Right-hand Side of the Kriging System", \
                             VAR, 2, "Variance of the Kriging System"

ENUM_DECLARE(ENUM_CALC_MEMBER)
