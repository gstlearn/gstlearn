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

class Db;
class DbGrid;

class GSTLEARN_EXPORT DbHelper
{
public:
  static int findDuplicates(Db *db1,
                            Db *db2,
                            bool flag_same,
                            bool verbose,
                            int opt_code,
                            double tolcode,
                            const VectorDouble &dist,
                            VectorDouble &sel);
  static int centerPointToGrid(Db *db_point, DbGrid *db_grid, double eps_random=EPSILON6);
  static int normalizeVariables(Db *db,
                                const char *oper,
                                const VectorInt& cols,
                                double center,
                                double stdv);
  static int dbgrid_filling(DbGrid *dbgrid,
                            int mode,
                            int seed,
                            int radius,
                            bool verbose = false,
                            const NamingConvention &namconv = NamingConvention("Fill"));
  static int db_duplicate(Db *db,
                          bool verbose = false,
                          const VectorDouble &dist = VectorDouble(),
                          int opt_code = 0,
                          double tolcode = 0.,
                          const NamingConvention &namconv = NamingConvention("Duplicate", true, true, true,
                                                                             ELoc::fromKey("SEL")));

  static int db_compositional_transform(Db *db,
                                        int verbose,
                                        int mode,
                                        int type,
                                        int number,
                                        int *iatt_in,
                                        int *iatt_out,
                                        int *numout);
  static DbGrid* dbgrid_sampling(DbGrid *dbin, const VectorInt &nmult);
  static int db_grid1D_fill(DbGrid *dbgrid,
                            int mode = 0,
                            int seed = 34243,
                            const NamingConvention &namconv = NamingConvention("Fill"));
};

//typedef VectorHelper VH;
class DbH: public DbHelper {};
