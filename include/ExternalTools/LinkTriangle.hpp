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

#include "triangle.hpp"
#include "Db/Db.hpp"
#include "Mesh/MeshEStandard.hpp"

GSTLEARN_EXPORT void meshes_1D_free(segmentio *t, int mode);
GSTLEARN_EXPORT void meshes_1D_init(segmentio *t);
GSTLEARN_EXPORT int meshes_1D_from_db(Db *db,segmentio *t);
GSTLEARN_EXPORT int meshes_1D_from_points(int nech, double *x, segmentio *t);
GSTLEARN_EXPORT void meshes_1D_default(Db *dbin, Db *dbout, segmentio *t);
GSTLEARN_EXPORT void meshes_1D_print(segmentio *t, int brief);
GSTLEARN_EXPORT void meshes_1D_create(int verbose,
                                      struct segmentio *in,
                                      struct segmentio *out);
GSTLEARN_EXPORT MeshEStandard* meshes_1D_load_vertices(segmentio *t);
GSTLEARN_EXPORT void meshes_1D_extended_domain(Db *dbout,
                                               const double *gext,
                                               segmentio *t);

GSTLEARN_EXPORT void meshes_2D_free(triangulateio *t, int mode);
GSTLEARN_EXPORT void meshes_2D_init(triangulateio *t);
GSTLEARN_EXPORT int meshes_2D_from_db(Db *db,
                                      int use_code,
                                      triangulateio *t);
GSTLEARN_EXPORT int meshes_2D_from_points(int nech,
                                          double *x,
                                          double *y,
                                          triangulateio *t);
GSTLEARN_EXPORT void meshes_2D_default(Db *dbin, Db *dbout, triangulateio *t);
GSTLEARN_EXPORT int meshes_2D_from_mem(int nseg,
                                       int ncol,
                                       int *segments,
                                       triangulateio *t);
GSTLEARN_EXPORT void meshes_2D_print(triangulateio *t, int brief);
GSTLEARN_EXPORT void meshes_2D_create(int verbose,
                                      const String& triswitches,
                                      struct triangulateio *in,
                                      struct triangulateio *out,
                                      struct triangulateio *vorout);
GSTLEARN_EXPORT MeshEStandard* meshes_2D_load_vertices(triangulateio *t);
GSTLEARN_EXPORT void meshes_2D_extended_domain(Db *dbout,
                                               const double *gext,
                                               triangulateio *t);
