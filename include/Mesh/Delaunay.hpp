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

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Mesh/AMesh.hpp"

namespace gstlrn
{
GSTLEARN_EXPORT VectorDouble get_db_extension(Db* dbin, Db* dbout, Id* nout);
GSTLEARN_EXPORT VectorDouble extend_grid(DbGrid* db, const VectorDouble& gext, Id* nout);
GSTLEARN_EXPORT VectorDouble extend_point(Db* db, const VectorDouble& gext, Id* nout);
GSTLEARN_EXPORT Id MSS(Id ndim, Id ipol, Id icas, Id icorn, Id idim);
GSTLEARN_EXPORT Id meshes_2D_write(const char* file_name,
                                    const char* obj_name,
                                    Id verbose,
                                    Id ndim,
                                    Id ncode,
                                    Id ntri,
                                    Id npoints,
                                    const VectorInt& ntcode,
                                    const VectorInt& triangles,
                                    const VectorDouble& points);
GSTLEARN_EXPORT AMesh* meshes_turbo_1D_grid_build(DbGrid* dbgrid);
GSTLEARN_EXPORT AMesh* meshes_turbo_2D_grid_build(DbGrid* dbgrid);
GSTLEARN_EXPORT AMesh* meshes_turbo_3D_grid_build(DbGrid* dbgrid);

GSTLEARN_EXPORT void mesh_stats(Id ndim,
                                Id ncorner,
                                Id nmesh,
                                const I32* meshes,
                                const double* points);
                                
}
