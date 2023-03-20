/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
