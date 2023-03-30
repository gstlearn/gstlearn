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

#define ENUM_ANAM EAnam, HERMITIAN, \
                  UNKNOWN,     0,  "Undefined", \
                  EXTERNAL,    1,  "External anamorphosis", \
                  HERMITIAN,   2,  "Hermitian anamorphosis", \
                  EMPIRICAL,   3,  "Empirical anamorphosis", \
                  DISCRETE_DD, 4,  "Disjunctive Discrete anamorphosis", \
                  DISCRETE_IR, 5,  "Indicator Residuals anamorphosis"

ENUM_DECLARE(ENUM_ANAM)
