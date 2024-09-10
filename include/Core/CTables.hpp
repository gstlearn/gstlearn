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
#include "geoslib_d.h"

GSTLEARN_EXPORT CTables* ct_tables_manage(int mode,
                                          int verbose,
                                          int flag_cumul,
                                          int nconf,
                                          int ndisc,
                                          double cmin,
                                          double cmax,
                                          CTables* ctables_old);
GSTLEARN_EXPORT void ct_tables_print(CTables* ctables, int flag_print);
GSTLEARN_EXPORT int ct_tableone_covrank(const CTables* ctables, double cova, double* cround);
GSTLEARN_EXPORT int ct_tableone_getrank_from_proba(CTables* ctables,
                                                   double gaussian);
GSTLEARN_EXPORT double ct_tableone_calculate(CTables* ctables, int iconf0, double* lows, double* ups);
GSTLEARN_EXPORT double ct_tableone_calculate_by_rank(CTables* ctables,
                                                     int iconf0,
                                                     double* rklows,
                                                     double* rkups);
GSTLEARN_EXPORT double ct_INTRES2(CTables* ctables, int iconf0, int idisc0, int jdisc0);
GSTLEARN_EXPORT double ct_INTRES3(CTables* ctables, int iconf0, int idisc0, int jdisc0, int kdisc0);
