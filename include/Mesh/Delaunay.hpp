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

GSTLEARN_EXPORT VectorDouble get_db_extension(Db* dbin, Db* dbout, int* nout);
GSTLEARN_EXPORT VectorDouble extend_grid(DbGrid* db, const VectorDouble& gext, int* nout);
GSTLEARN_EXPORT VectorDouble extend_point(Db* db, const VectorDouble& gext, int* nout);
GSTLEARN_EXPORT int MSS(int ndim, int ipol, int icas, int icorn, int idim);
GSTLEARN_EXPORT int meshes_2D_write(const char* file_name,
                                    const char* obj_name,
                                    int verbose,
                                    int ndim,
                                    int ncode,
                                    int ntri,
                                    int npoints,
                                    const VectorInt& ntcode,
                                    const VectorInt& triangles,
                                    const VectorDouble& points);
GSTLEARN_EXPORT AMesh* meshes_turbo_1D_grid_build(DbGrid* dbgrid);
GSTLEARN_EXPORT AMesh* meshes_turbo_2D_grid_build(DbGrid* dbgrid);
GSTLEARN_EXPORT AMesh* meshes_turbo_3D_grid_build(DbGrid* dbgrid);

GSTLEARN_EXPORT void mesh_stats(int ndim,
                                int ncorner,
                                int nmesh,
                                const int* meshes,
                                const double* points);
                                
