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

#define ENUM_SELECTIVITY ESelectivity, UNKNOWN, \
                 UNKNOWN, -1, "Unknown Option", \
                 Z,        0, "Grade", \
                 T,        1, "Tonnage", \
                 Q,        2, "Metal quantity", \
                 B,        3, "Conventional Benefit", \
                 M,        4, "Recovered mean", \
                 PROP,     5, "Probability to exceed Cutoff", \
                 QUANT,    6, "Quantile"

ENUM_DECLARE(ENUM_SELECTIVITY)

