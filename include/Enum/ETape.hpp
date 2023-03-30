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

#define ENUM_TAPE ETape, SPHERICAL,\
                 SPHERICAL,   0, "Spherical",\
                 CUBIC,       1, "Cubic",\
                 TRIANGLE,    2, "Triangle",\
                 PENTAMODEL,  3, "PentaModel",\
                 STORKEY,     4, "Storkey",\
                 WENDLAND1,   5, "Wendland-1",\
                 WENDLAND2,   6, "Wendland-2"

ENUM_DECLARE(ENUM_TAPE)
