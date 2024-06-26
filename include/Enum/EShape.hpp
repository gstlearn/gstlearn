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

#define ENUM_SHAPE EShape, PARALLELEPIPED,\
                 PARALLELEPIPED,    0, "Parallelepiped",\
                 ELLIPSOID,         1, "Full Ellipsoid",\
                 PARABOLOID,        2, "Full Paraboloid",\
                 HALFELLIPSOID,     3, "Lower-Half Ellipsoid",\
                 HALFPARABOLOID,    4, "Lower-Hald Paraboloid",\
                 HALFSINUSOID,      5, "Lower-Half Sinusoid"

ENUM_DECLARE(ENUM_SHAPE)

