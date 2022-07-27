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

#define ENUM_SELECTIVITY ESelectivity, UNKNOWN, \
                 UNKNOWN, -1, "Unknown Option", \
                 Z,        0, "Grade", \
                 T,        1, "Tonnage", \
                 Q,        2, "Metal quantity", \
                 B,        3, "Conventional Benefit", \
                 M,        4, "Recovered mean", \
                 PROBA,    5, "Probability to exceed Cutoff", \
                 QUANT,    6, "Quantile"

ENUM_DECLARE(ENUM_SELECTIVITY)

