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

#define ENUM_CONVDIR EConvDir, X,\
                 X,           0, "Along X",\
                 Y,           1, "Along Y",\
                 Z,           2, "Along Z",\
                 XY,          3, "Along X-Y",\
                 XYZ,         4, "Along X-Y-Z"

ENUM_DECLARE(ENUM_CONVDIR)
