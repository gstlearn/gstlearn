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

#define ENUM_TSHAPE ETShape, PARALLELEPIPED,\
                 PARALLELEPIPED,    0, "Parallelepiped",\
                 LOWHALFELLIPSOID,  1, "Lower-Half Ellipsoid",\
                 UPHALFELLIPSOID,   2, "Upper-Half Ellipsoid",\
                 LOWHALFPARABOLOID, 3, "Lower-Hald Paraboloid",\
                 UPHALFPARABOLOID,  4, "Upper-Half Paraboloid",\
                 UPHALFSINUSOID,    5, "Upper-Half Sinusoid",\
                 ELLPSOID,          6, "Full Ellipsoid",\
                 PARABOLOID,        7, "Full Paraboloid"

ENUM_DECLARE(ENUM_TSHAPE)

