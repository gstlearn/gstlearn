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

#include "Basic/NamingConvention.hpp"
#include "geoslib_define.h"

#include "Enum/ECalcVario.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class Vario;
  class Db;
  class DbGrid;
  class VarioParam;
  class Polygons;

  GSTLEARN_EXPORT Vario* variogramCalculate(
    Db* db,
    const ECalcVario& calculType = ECalcVario::fromKey("VARIOGRAM"),
    bool flag_ergodic = true,
    Id nlag = 10,
    double dlag = 1.,
    Id ndir = 1,
    const VectorDouble& angles = VectorDouble(),
    double toldis = 0.5,
    double tolang = TEST,
    bool verbose = false);

  GSTLEARN_EXPORT Vario* varioGridCalculate(
    DbGrid* dbgrid,
    const ECalcVario& calculType = ECalcVario::fromKey("VARIOGRAM"),
    bool flag_ergodic = true,
    Id nlag = 10,
    bool flagAllDirections = true,
    const VectorVectorInt& dirincr = VectorVectorInt(),
    bool verbose = false);

  GSTLEARN_EXPORT VectorInt variogramPerPoint(
    Db* db,
    Id iech0,
    Id ilag0 = 10,
    double dlag = 1.,
    const VectorDouble& codir = VectorDouble(),
    double toldis = 0.5,
    double tolang = TEST,
    double bench = TEST,
    double cylrad = TEST,
    const VectorDouble& benchdir = VectorDouble());

  GSTLEARN_EXPORT DbGrid* vmapFromDb(
    Db* db,
    Id radius = 0,
    bool flag_FFT = true,
    const ECalcVario& calculType = ECalcVario::fromKey("VARIOGRAM"),
    bool flag_ergodic = true,
    const VectorInt& nxx = VectorInt(),
    const VectorDouble& dxx = VectorDouble(),
    const NamingConvention& namconv = NamingConvention("VMAP"));

  GSTLEARN_EXPORT DbGrid* vcloudFromDb(
    Db* db,
    const VarioParam* varioparam = nullptr,
    double lagmax = TEST,
    const ECalcVario& calculType = ECalcVario::fromKey("VARIOGRAM"),
    bool flag_ergodic = true,
    double varmax = TEST,
    Id lagnb = 100,
    Id varnb = 100,
    const Polygons* polygon = nullptr,
    double distmax = TEST,
    double varmin = TEST,
    Id countmax = ITEST,
    const VectorInt& samples = VectorInt(),
    const NamingConvention& namconv = NamingConvention("Cloud"));

  GSTLEARN_EXPORT DbGrid* vcloudFromParam(
    Db* db,
    const ECalcVario& calculType = ECalcVario::fromKey("VARIOGRAM"),
    bool flag_ergodic = true,
    Id nlag = 10,
    double dlag = 1.,
    Id ndir = 1,
    const VectorDouble& angles = VectorDouble(),
    double toldis = 0.5,
    double tolang = TEST,
    Id lagnb = 100,
    Id varnb = 100,
    const Polygons* polygon = nullptr,
    double distmax = TEST,
    double varmin = TEST,
    Id countmax = ITEST,
    const VectorInt& samples = VectorInt());

} // namespace gstlrn
