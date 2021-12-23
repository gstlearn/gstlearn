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

#define ENUM_TAPE ETape, SPHERICAL,\
                 SPHERICAL,   0, "Spherical",\
                 CUBIC,       1, "Cubic",\
                 TRIANGLE,    2, "Triangle",\
                 PENTAMODEL,  3, "PentaModel",\
                 STORKEY,     4, "Storkey",\
                 WENDLAND1,   5, "Wendland-1",\
                 WENDLAND2,   6, "Wendland-2"

ENUM_DECLARE(ENUM_TAPE)
