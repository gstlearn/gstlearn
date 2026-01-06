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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Variogram/AVario.hpp"
#include "Variogram/VarioParam.hpp"

namespace gstlrn
{
class Db;
class ECalcVario;
class Polygons;

/**
 * \brief
 * Class containing the Variogram Cloud which uses an DbGrid provided by the user
 * This function simply calculate and add the results as new field in this DbGrid.
 */
class GSTLEARN_EXPORT VCloud: public AVario
{
public:
  VCloud(DbGrid* dbcloud, const VarioParam* varioparam);
  VCloud(const VCloud& r);
  VCloud& operator=(const VCloud& r);
  virtual ~VCloud();

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(VCloud)

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AVCloud Interface
  double _getIVAR(const Db* db, Id iech, Id ivar) const override;
  void _setAVarioResult(Id iech1,
                        Id iech2,
                        Id nvar,
                        Id idir,
                        Id ilag,
                        Id ivar,
                        Id jvar,
                        Id orient,
                        double ww,
                        double w1,
                        double w2,
                        double z1,
                        double z2,
                        double dist,
                        double value) override;

  Id compute(Db* db,
             const ECalcVario& calculType    = ECalcVario::fromKey("VARIOGRAM"),
             bool flag_ergodic               = true,
             const NamingConvention& namconv = NamingConvention("Cloud"));

  Id selectFromPolygon(Db* db, Polygons* polygon, Id idir = 0);

  void dumpStorage(Id mode        = 0,
                   double distmax = TEST,
                   double varmin  = 0.,
                   Id npairmax    = ITEST,
                   Id ndataMax    = ITEST);

private:
  void _variogram_cloud(Db* db, Id idir);
  void _final_discretization_grid();
  Id _update_discretization_grid(double x, double y);

private:
  DbGrid* _dbcloud;              // Pointer to the already existing output DbGrid (not to be deleted)
  const VarioParam* _varioparam; // Pointer (not to be deleted)
  Id _IPTR;
};

GSTLEARN_EXPORT DbGrid* vcloudGrid(const Db* db,
                                   double lagmax = TEST,
                                   double varmax = TEST,
                                   Id lagnb      = 100,
                                   Id varnb      = 100);

GSTLEARN_EXPORT DbGrid* db_vcloud(Db* db,
                                  const VarioParam* varioparam,
                                  const ECalcVario& calculType    = ECalcVario::fromKey("VARIOGRAM"),
                                  bool flag_ergodic               = true,
                                  double lagmax                   = TEST,
                                  double varmax                   = TEST,
                                  Id lagnb                        = 100,
                                  Id varnb                        = 100,
                                  const NamingConvention& namconv = NamingConvention("Cloud"));

GSTLEARN_EXPORT DbGrid* vcloudCalculate(Db* db,
                                        const ECalcVario& calculType = ECalcVario::fromKey("VARIOGRAM"),
                                        bool flag_ergodic            = true,
                                        Id nlag                      = 10,
                                        double dlag                  = 1.,
                                        Id ndir                      = 1,
                                        const VectorDouble& angles   = VectorDouble(),
                                        double toldis                = 0.5,
                                        double tolang                = TEST,
                                        Id lagnb                     = 100,
                                        Id varnb                     = 100);
} // namespace gstlrn