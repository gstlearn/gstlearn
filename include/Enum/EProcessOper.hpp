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

#define ENUM_PROCESS_OPER EProcessOper, UNDEFINED, \
                          UNDEFINED,  -1, "Undefined Process Operation", \
                          COPY,        0, "Copy Process", \
                          MARGINAL,    1, "Marginal Process", \
                          CONDITIONAL, 2, "Conditional Process"

ENUM_DECLARE(ENUM_PROCESS_OPER)
