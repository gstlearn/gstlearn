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

#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

class CSVformat;
class Db;

GSTLEARN_EXPORT int csv_manage(const char* filename,
                               const CSVformat& csv,
                               int mode,
                               int nitem,
                               bool flagInteger = false,
                               bool verbose     = false);
GSTLEARN_EXPORT void csv_print_double(double value);

GSTLEARN_EXPORT Db* db_read_csv(const char* file_name,
                                const CSVformat& csvfmt,
                                int verbose            = 0,
                                int ncol_max           = -1,
                                int nrow_max           = -1,
                                bool flagAddSampleRank = false);
GSTLEARN_EXPORT int db_write_csv(Db* db,
                                 const char* filename,
                                 const CSVformat& csv,
                                 int flag_allcol  = 1,
                                 int flag_coor    = 1,
                                 bool flagInteger = false);
GSTLEARN_EXPORT int csv_table_read(const String& filename,
                                   const CSVformat& csvfmt,
                                   int verbose,
                                   int ncol_max,
                                   int nrow_max,
                                   int* ncol_arg,
                                   int* nrow_arg,
                                   VectorString& names,
                                   VectorDouble& tab);
