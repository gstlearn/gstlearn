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
                 NUM,      0, "Number of defined values", \
                 MEAN,     1, "Mean over the defined values", \
                 VAR,      2, "Variance over the defined values", \
                 STDV,     3, "Standard Deviation over the defined values", \
                 MINI,     4, "Minimum over the defined values", \
                 MAXI,     5, "Maximum over the defined values", \
                 SUM,      6, "Sum over the variables", \
                 PROP,     7, "Proportion of values within [vmin;vmax]", \
                 QUANT,    8, "Quantile corresponding to given probability", \
                 T,        9, "Tonnage within [vmin;vmax]", \
                 Q,       10, "Metal quantity within [vmin;vmax]", \
                 M,       11, "Recovered mean within [vmin;vmax]", \
                 B,       12, "Conventional Benefit within [vmin;vmax]", \
                 COV,     13, "Covariance", \
                 CORR,    14, "Correlation", \
                 ZERO,    15, "Count of zero values", \
                 MEDIAN,  16, "Median", \
                 MEAN2,   17, "Mean of the defined values for secondary variable", \
                 VAR2,    18, "Variance over the defined values f secondary variable", \
                 STDV2,   19, "Standard Deviation for the secondary variable", \
                 SUM2,    20, "Sum over the secondary variable", \
                 PLUS,    21, "Count of positive values", \
                 MOINS,   22, "Count of negative values", \
                 ORE,     23, "Ore quantity", \
                 METAL,   24, "Metal quantity"


ENUM_DECLARE(ENUM_STATOPTION)
