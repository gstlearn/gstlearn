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

#define ENUM_ANAM EAnam, HERMITIAN, \
                  UNKNOWN,     0,  "Undefined", \
                  EXTERNAL,    1,  "External anamorphosis", \
                  HERMITIAN,   2,  "Hermitian anamorphosis", \
                  EMPIRICAL,   3,  "Empirical anamorphosis", \
                  DISCRETE_DD, 4,  "Disjunctive Discrete anamorphosis", \
                  DISCRETE_IR, 5,  "Indicator Residuals anamorphosis"

ENUM_DECLARE(ENUM_ANAM)
