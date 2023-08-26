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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

class SphTriangle;
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
