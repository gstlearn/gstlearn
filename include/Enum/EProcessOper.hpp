/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_PROCESS_OPER EProcessOper, UNDEFINED, \
                          UNDEFINED,  -1, "Undefined Process Operation", \
                          COPY,        0, "Copy Process", \
                          MARGINAL,    1, "Marginal Process", \
                          CONDITIONAL, 2, "Conditional Process"

ENUM_DECLARE(ENUM_PROCESS_OPER)
