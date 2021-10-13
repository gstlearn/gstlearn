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

#define ENUM_ANAM EAnam, UNDEFINED, \
                  UNDEFINED,  -1,  "Undefined anamorphosis", \
                  EXTERNAL,    0,  "External anamorphosis", \
                  HERMITIAN,   1,  "Hermitian anamorphosis", \
                  EMPIRICAL,   2,  "Empirical anamorphosis", \
                  DISCRETE_DD, 3,  "Discrete anamorphosis", \
                  DISCRETE_IR, 4,  "Discrete Indicator Residuals anamorphosis"

ENUM_DECLARE(ENUM_ANAM)
