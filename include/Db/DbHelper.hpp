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

#include "Basic/VectorNumT.hpp"
#include "Basic/NamingConvention.hpp"

namespace gstlrn
{
class Db;
class DbGrid;

class GSTLEARN_EXPORT DbHelper
{
public:
  static Id findDuplicates(Db *db1,
                            Db *db2,
                            bool flag_same,
                            bool verbose,
                            Id opt_code,
                            double tolcode,
                            const VectorDouble &dist,
                            VectorDouble &sel);
  static Id centerPointToGrid(Db *db_point, DbGrid *db_grid, double eps_random=EPSILON6);
  static Id normalizeVariables(Db *db,
                                const char *oper,
                                const VectorInt& cols,
                                double center,
                                double stdv);
  static Id dbgrid_filling(DbGrid *dbgrid,
                            Id mode = 0,
                            Id seed = 34342,
                            Id radius = 1,
                            bool verbose = false,
                            const NamingConvention &namconv = NamingConvention("Fill"));
  static Id db_duplicate(Db *db,
                          bool verbose = false,
                          const VectorDouble &dist = VectorDouble(),
                          Id opt_code = 0,
                          double tolcode = 0.,
                          const NamingConvention &namconv = NamingConvention("Duplicate", true, true, true,
                                                                             ELoc::fromKey("SEL")));

  static Id db_compositional_transform(Db *db,
                                        Id verbose,
                                        Id mode,
                                        Id type,
                                        Id number,
                                        Id *iatt_in,
                                        Id *iatt_out,
                                        Id *numout);
  static DbGrid* dbgrid_sampling(DbGrid *dbin, const VectorInt &nmult);
  static Id db_grid1D_fill(DbGrid *dbgrid,
                            Id mode = 0,
                            Id seed = 34243,
                            const NamingConvention &namconv = NamingConvention("Fill"));
};

//typedef VectorHelper VH;
class DbH: public DbHelper {};
}