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

namespace gstlrn
{
class CSVformat;
class Db;

GSTLEARN_EXPORT Id csv_manage(const char* filename,
                               const CSVformat& csv,
                               Id mode,
                               Id nitem,
                               bool flagInteger = false,
                               bool verbose     = false);
GSTLEARN_EXPORT void csv_print_double(double value);

GSTLEARN_EXPORT Db* db_read_csv(const String& filename,
                                const CSVformat& csvfmt,
                                bool verbose           = false,
                                Id ncol_max           = -1,
                                Id nrow_max           = -1,
                                bool flagAddSampleRank = false);
GSTLEARN_EXPORT Id db_write_csv(Db* db,
                                 const char* filename,
                                 const CSVformat& csv,
                                 Id flag_allcol  = 1,
                                 Id flag_coor    = 1,
                                 bool flagInteger = false);
GSTLEARN_EXPORT Id csv_table_read(const String& filename,
                                   const CSVformat& csvfmt,
                                   bool verbose,
                                   Id ncol_max,
                                   Id nrow_max,
                                   Id* ncol_arg,
                                   Id* nrow_arg,
                                   VectorString& names,
                                   VectorDouble& tab);
} // namespace gstlrn