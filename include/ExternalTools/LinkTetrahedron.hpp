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

#include "tetgen.hpp"

class Db;
class MeshEStandard;

GSTLEARN_EXPORT void meshes_3D_free(tetgenio *t);
GSTLEARN_EXPORT void meshes_3D_create(int verbose,
                                      const String &triswitch,
                                      tetgenio *in,
                                      tetgenio *out);
GSTLEARN_EXPORT void meshes_3D_print(tetgenio *t, int brief);
GSTLEARN_EXPORT int meshes_3D_from_db(Db *db,tetgenio *t);
GSTLEARN_EXPORT int meshes_3D_from_points(int nech,
                                          double *x,
                                          double *y,
                                          double *z,
                                          tetgenio *t);
GSTLEARN_EXPORT void meshes_3D_default(Db *dbin, Db *dbout, tetgenio *t);
GSTLEARN_EXPORT MeshEStandard* meshes_3D_load_vertices(tetgenio *t);
GSTLEARN_EXPORT void meshes_3D_extended_domain(Db *dbout,
                                               const double *gext,
                                               tetgenio *t);
