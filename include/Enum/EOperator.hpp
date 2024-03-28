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

#define ENUM_OPERATOR EOperator, IDLE, \
                 IDLE,      0, "New = New",\
                 ADD,       1, "New = New + Old", \
                 PRODUCT,   2, "New = New * Old", \
                 SUBTRACT,  3, "New = New - Old", \
                 SUBOPP,    4, "New = Old - New", \
                 DIVIDE,    5, "New = New / Old", \
                 DIVOPP,    6, "New = Old / New", \
                 DEFINE,    7, "New = New (if Old is defined) or TEST otherwise ", \
                 MIN,       8, "New = MIN(New, Old)", \
                 MAX,       9, "New = MAX(New, Old)"

ENUM_DECLARE(ENUM_OPERATOR)
