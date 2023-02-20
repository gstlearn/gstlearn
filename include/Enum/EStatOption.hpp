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

#define ENUM_STATOPTION EStatOption, UNKNOWN, \
                 UNKNOWN, -1, "Unknown Option", \
                 NUM,      0, "Number", \
                 MEAN,     1, "Mean", \
                 VAR,      2, "Variance", \
                 STDV,     3, "St. Dev.", \
                 MINI,     4, "Minimum", \
                 MAXI,     5, "Maximum", \
                 SUM,      6, "Sum", \
                 PROP,     7, "Prop.", \
                 QUANT,    8, "Quantile", \
                 T,        9, "Tonnage", \
                 Q,       10, "Metal", \
                 M,       11, "Rec. Mean", \
                 B,       12, "Benefit", \
                 COV,     13, "Covariance", \
                 CORR,    14, "Correlation", \
                 ZERO,    15, "Zero Count", \
                 MEDIAN,  16, "Median", \
                 MEAN2,   17, "Mean of the defined values for secondary variable", \
                 VAR2,    18, "Variance over the defined values f secondary variable", \
                 STDV2,   19, "Standard Deviation for the secondary variable", \
                 SUM2,    20, "Sum over the secondary variable", \
                 PLUS,    21, "Count of positive values", \
                 MOINS,   22, "Count of negative values", \
                 ORE,     23, "Ore", \
                 METAL,   24, "Metal"


ENUM_DECLARE(ENUM_STATOPTION)
