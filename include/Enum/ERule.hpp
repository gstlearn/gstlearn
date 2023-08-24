/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_RULE ERule, STD, \
                  STD,    0,  "Standard Lithotype Rule", \
                  SHIFT,  1,  "Shift Lithotype Rule", \
                  SHADOW, 2,  "Shadow Lithotype Rule"

ENUM_DECLARE(ENUM_RULE)
