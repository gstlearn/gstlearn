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

class DbGrid;
class Model;

GSTLEARN_EXPORT int seismic_estimate_XZ(DbGrid* db,
                                        Model* model,
                                        int nbench,
                                        int nv2max,
                                        int flag_ks,
                                        int flag_std,
                                        int flag_sort,
                                        int flag_stat);
GSTLEARN_EXPORT int seismic_simulate_XZ(DbGrid* db,
                                        Model* model,
                                        int nbench,
                                        int nv2max,
                                        int nbsimu,
                                        int seed,
                                        int flag_ks,
                                        int flag_sort,
                                        int flag_stat);
GSTLEARN_EXPORT int seismic_z2t_grid(
  int verbose, DbGrid* db_z, int iatt_v, int* nx, double* x0, double* dx);
GSTLEARN_EXPORT int seismic_t2z_grid(
  int verbose, DbGrid* db_t, int iatt_v, int* nx, double* x0, double* dx);
GSTLEARN_EXPORT int seismic_z2t_convert(DbGrid* db_z, int iatt_v, DbGrid* db_t);
GSTLEARN_EXPORT int seismic_t2z_convert(DbGrid* db_t, int iatt_v, DbGrid* db_z);
GSTLEARN_EXPORT int seismic_operate(DbGrid* db, int oper);
GSTLEARN_EXPORT int seismic_convolve(DbGrid* db,
                                     int flag_operate,
                                     int flag_contrast,
                                     int type,
                                     int ntw,
                                     int option,
                                     int tindex,
                                     double fpeak,
                                     double period,
                                     double amplitude,
                                     double distort,
                                     double val_before,
                                     double val_middle,
                                     double val_after,
                                     double* wavelet);
