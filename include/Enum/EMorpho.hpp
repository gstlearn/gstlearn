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

#define ENUM_MORPHO EMorpho, UNKNOWN,\
                     UNKNOWN,   0, "Idle", \
                     THRESH,    1, "Convert the Input Variable into Binary Image", \
                     NEGATION,  2, "Invert of the Binary Image", \
                     EROSION,   3, "Erosion on the Binary Image", \
                     DILATION,  4, "Dilation on the Binary Image", \
                     OPEN,      5, "Opening (erosion then dilation) on the Binary Image", \
                     CLOSE,     6, "Closing (dilation then erosion) on the Binary Image", \
                     CC,        7, "Connected components (cells assigned Rank of CC)", \
                     CCSIZE,    8, "Connected components (cells assigned Volume of CC)", \
                     DISTANCE,  9, "Distance to the pore edge", \
                     ANGLE,    10, "Angle of the tangent to isovalues of a coloured image", \
                     GRADIENT, 11, "Gradient components"

ENUM_DECLARE(ENUM_MORPHO)
