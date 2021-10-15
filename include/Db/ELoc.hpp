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

#define ENUM_LOC ELoc, UNKNOWN, \
                 UNKNOWN, -1, "Unknown locator", \
                 X,        0, "Coordinate", \
                 Z,        1, "Variable", \
                 V,        2, "Variance of measurement error", \
                 F,        3, "External Drift", \
                 G,        4, "Gradient component", \
                 L,        5, "Lower bound of an inequality", \
                 U,        6, "Upper bound of an inequality", \
                 P,        7, "Proportion", \
                 W,        8, "Weight", \
                 C,        9, "Code", \
                 SEL,     10, "Selection", \
                 DOM,     11, "Domain", \
                 BLEX,    12, "Block Extension", \
                 ADIR,    13, "Dip direction Angle", \
                 ADIP,    14, "Dip Angle", \
                 SIZE,    15, "Object height", \
                 BU,      16, "Fault UP termination", \
                 BD,      17, "Fault DOWN termination", \
                 TIME,    18, "Time variable", \
                 LAYER,   19, "Layer rank", \
                 NOSTAT,  20, "Non-stationary parameter", \
                 TGTE,    21, "Tangent", \
                 SIMU,    22, "Conditional or non-conditional simulations", \
                 FACIES,  23, "Facies simulated", \
                 GAUSFAC, 24, "Gaussian value for Facies", \
                 DATE,    25, "Date", \
                 RKLOW,   26, "Rank for lower bound (when discretized)", \
                 RKUP,    27, "Rank for upper bound (when discretized)"

ENUM_DECLARE(ENUM_LOC)
