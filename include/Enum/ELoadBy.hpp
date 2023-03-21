/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_LOAD_BY ELoadBy, COLUMN,\
                     COLUMN,   0, "Values are provided sorted by columns", \
                     SAMPLE,   1, "Values are provided sorted by sample"

ENUM_DECLARE(ENUM_LOAD_BY)
