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

#include "AEnum.hpp"

#define ENUM_NEIGH ENeigh, UNIQUE, \
                   UNIQUE, 0, "Unique Neighborhood", \
                   BENCH,  1, "Bench Neighborhood", \
                   MOVING, 2, "Moving Neighborhood", \
                   IMAGE,  3, "Image Neighborhood"

ENUM_DECLARE(ENUM_NEIGH)
