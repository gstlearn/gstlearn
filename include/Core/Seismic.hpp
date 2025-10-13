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
#include "gstlearn_export.hpp"

namespace gstlrn
{
class DbGrid;
class Model;

GSTLEARN_EXPORT Id seismic_estimate_XZ(DbGrid* db,
                                        Model* model,
                                        Id nbench,
                                        Id nv2max,
                                        Id flag_ks,
                                        Id flag_std,
                                        Id flag_sort,
                                        Id flag_stat);
GSTLEARN_EXPORT Id seismic_simulate_XZ(DbGrid* db,
                                        Model* model,
                                        Id nbench,
                                        Id nv2max,
                                        Id nbsimu,
                                        Id seed,
                                        Id flag_ks,
                                        Id flag_sort,
                                        Id flag_stat);
GSTLEARN_EXPORT Id seismic_z2t_grid(
  Id verbose,
  DbGrid* db_z,
  Id iatt_v,
  Id* nx,
  double* x0,
  double* dx);
GSTLEARN_EXPORT Id seismic_t2z_grid(
  Id verbose,
  DbGrid* db_t,
  Id iatt_v,
  Id* nx,
  double* x0,
  double* dx);
GSTLEARN_EXPORT Id seismic_z2t_convert(DbGrid* db_z, Id iatt_v, DbGrid* db_t);
GSTLEARN_EXPORT Id seismic_t2z_convert(DbGrid* db_t, Id iatt_v, DbGrid* db_z);
GSTLEARN_EXPORT Id seismic_operate(DbGrid* db, Id oper);
GSTLEARN_EXPORT Id seismic_convolve(DbGrid* db,
                                     Id flag_operate,
                                     Id flag_contrast,
                                     Id type,
                                     Id ntw,
                                     Id option,
                                     Id tindex,
                                     double fpeak,
                                     double period,
                                     double amplitude,
                                     double distort,
                                     double val_before,
                                     double val_middle,
                                     double val_after,
                                     VectorDouble& wavelet);
}
