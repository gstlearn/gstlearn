/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
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

