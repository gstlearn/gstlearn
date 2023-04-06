/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_JUSTIFY EJustify, LEFT, \
                     LEFT,  -1,  "Left Justification", \
                     CENTER, 0,  "Center Justification", \
                     RIGHT,  1,  "Right Justification"

ENUM_DECLARE(ENUM_JUSTIFY)
