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

#define ENUM_POST_STAT EPostStat, UNKNOWN, \
                 UNKNOWN, -1, "Unknown Option", \
                 MEAN,     1, "Mean", \
                 VAR,      2, "Variance", \
                 VARP,     3, "Variance-P", \
                 STD,      4, "Std.", \
                 STDP,     5, "Std-P", \
                 MED,      6, "Median", \
                 MINI,     7, "Minimum", \
                 MAXI,     8, "Maximum"

ENUM_DECLARE(ENUM_POST_STAT)
