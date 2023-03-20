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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include "Mesh/sphtriangle.hpp"

class Db;

GSTLEARN_EXPORT void meshes_2D_sph_init(SphTriangle *t);
GSTLEARN_EXPORT void meshes_2D_sph_free(SphTriangle *t, int mode);
GSTLEARN_EXPORT int meshes_2D_sph_from_db(Db *db, SphTriangle *t);
GSTLEARN_EXPORT int meshes_2D_sph_from_points(int nech,
                                              double *x,
                                              double *y,
                                              SphTriangle *t);
GSTLEARN_EXPORT int meshes_2D_sph_from_auxiliary(const String &triswitch,
                                                 SphTriangle *t);
GSTLEARN_EXPORT void meshes_2D_sph_print(SphTriangle *t, int brief);
GSTLEARN_EXPORT int meshes_2D_sph_create(int verbose, SphTriangle *t);
