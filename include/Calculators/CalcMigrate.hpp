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

#include "ACalcDbToDb.hpp"

#include "geoslib_define.h"

namespace gstlrn
{
class Db;
class DbGrid;
class EStatOption;

class GSTLEARN_EXPORT CalcMigrate: public ACalcDbToDb
{
public:
  CalcMigrate();
  CalcMigrate(const CalcMigrate& r)            = delete;
  CalcMigrate& operator=(const CalcMigrate& r) = delete;
  virtual ~CalcMigrate();

  void setDmax(const VectorDouble& dmax) { _dmax = dmax; }
  void setFlagFill(bool flagFill) { _flagFill = flagFill; }
  void setFlagInter(bool flagInter) { _flagInter = flagInter; }
  void setFlagBall(bool flagBall) { _flagBall = flagBall; }
  void setDistType(Id dist_type) { _distType = dist_type; }
  void setIuids(const VectorInt& iuids) { _iuids = iuids; }
  void setFlagLocate(bool flagLocate) { _flagLocate = flagLocate; }
  void setLocatorType(const ELoc& locatorType) { _locatorType = locatorType; }

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

  static Id _migrate(Db* db1,
                      Db* db2,
                      Id iatt1,
                      Id iatt2,
                      Id distType,
                      const VectorDouble& dmax,
                      bool flag_fill,
                      bool flag_inter,
                      bool flag_ball);
  static Id _migratePointToGrid(Db* db_point,
                                 DbGrid* db_grid,
                                 Id iatt,
                                 Id distType,
                                 const VectorDouble& dmax,
                                 VectorDouble& tab);
  static Id _expandPointToPointBall(Db* db1,
                                     Db* db2,
                                     Id iatt,
                                     Id distType,
                                     const VectorDouble& dmax,
                                     VectorDouble& tab);
  static Id _migrateGridToGrid(DbGrid* db_gridin,
                                DbGrid* db_gridout,
                                Id iatt,
                                Id distType,
                                const VectorDouble& dmax,
                                VectorDouble& tab);
  static Id _expandPointToPoint(Db* db1,
                                 Db* db2,
                                 Id iatt,
                                 Id distType,
                                 const VectorDouble& dmax,
                                 VectorDouble& tab);
  static Id _expandGridToGrid(DbGrid* db_gridin,
                               DbGrid* db_gridout,
                               Id iatt,
                               Id distType,
                               const VectorDouble& dmax,
                               VectorDouble& tab);
  static Id _interpolateGridToPoint(DbGrid* db_grid,
                                     Db* db_point,
                                     Id iatt,
                                     Id distType,
                                     const VectorDouble& dmax,
                                     VectorDouble& tab);
  static Id _migrateGridToPoint(DbGrid* db_grid,
                                 Db* db_point,
                                 Id iatt,
                                 Id distType,
                                 const VectorDouble& dmax,
                                 VectorDouble& tab);

private:
  Id _iattOut;
  VectorInt _iuids;
  Id _distType;
  VectorDouble _dmax;
  bool _flagFill;
  bool _flagInter;
  bool _flagLocate;
  bool _flagBall;
  ELoc _locatorType;
};

GSTLEARN_EXPORT Id migrate(Db* dbin,
                            Db* dbout,
                            const String& name,
                            Id dist_type                   = 1,
                            const VectorDouble& dmax        = VectorDouble(),
                            bool flag_fill                  = false,
                            bool flag_inter                 = false,
                            bool flag_ball                  = false,
                            const NamingConvention& namconv = NamingConvention(
                              "Migrate",
                              false));
GSTLEARN_EXPORT Id migrateMulti(Db* dbin,
                                 Db* dbout,
                                 const VectorString& names,
                                 Id dist_type                   = 1,
                                 const VectorDouble& dmax        = VectorDouble(),
                                 bool flag_fill                  = false,
                                 bool flag_inter                 = false,
                                 bool flag_ball                  = false,
                                 const NamingConvention& namconv = NamingConvention(
                                   "Migrate"));
GSTLEARN_EXPORT Id migrateByAttribute(Db* dbin,
                                       Db* dbout,
                                       const VectorInt& iatts          = VectorInt(),
                                       Id dist_type                   = 1,
                                       const VectorDouble& dmax        = VectorDouble(),
                                       bool flag_fill                  = false,
                                       bool flag_inter                 = false,
                                       bool flag_ball                  = false,
                                       const NamingConvention& namconv = NamingConvention(
                                         "Migrate"));
GSTLEARN_EXPORT Id migrateByLocator(Db* dbin,
                                     Db* dbout,
                                     const ELoc& locatorType,
                                     Id dist_type                   = 1,
                                     const VectorDouble& dmax        = VectorDouble(),
                                     bool flag_fill                  = false,
                                     bool flag_inter                 = false,
                                     bool flag_ball                  = false,
                                     const NamingConvention& namconv = NamingConvention(
                                       "Migrate"));
GSTLEARN_EXPORT Id manageExternalInformation(Id mode,
                                              const ELoc& locatorType,
                                              Db* dbin,
                                              Db* dbout,
                                              bool* flag_created);
GSTLEARN_EXPORT Id interpolateVariableToPoint(DbGrid* db_grid,
                                               Id iatt,
                                               Id np,
                                               const double* xp,
                                               const double* yp,
                                               const double* zp,
                                               double* tab);
GSTLEARN_EXPORT VectorDouble dbgridLineSampling(DbGrid* dbgrid,
                                                const double* x1,
                                                const double* x2,
                                                Id ndisc,
                                                Id ncut,
                                                const double* cuts,
                                                Id* nval_ret);
GSTLEARN_EXPORT Id expandPointToGrid(Db* db_point,
                                      DbGrid* db_grid,
                                      Id iatt,
                                      Id iatt_time,
                                      Id iatt_angle,
                                      Id iatt_scaleu,
                                      Id iatt_scalev,
                                      Id iatt_scalew,
                                      Id flag_index,
                                      Id distType,
                                      const VectorDouble& dmax,
                                      VectorDouble& tab);
GSTLEARN_EXPORT Id pointToBlock(Db* dbpoint,
                                 DbGrid* dbgrid,
                                 Id option,
                                 Id flag_size,
                                 Id iatt_time,
                                 Id iatt_size,
                                 Id iatt_angle,
                                 Id iatt_scaleu,
                                 Id iatt_scalev,
                                 Id iatt_scalew);
GSTLEARN_EXPORT Id migrateGridToCoor(const DbGrid* db_grid,
                                      Id iatt,
                                      const VectorVectorDouble& coords,
                                      VectorDouble& tab);
GSTLEARN_EXPORT Id expandPointToCoor(const Db* db1,
                                      Id iatt,
                                      const VectorVectorDouble& coords,
                                      VectorDouble& tab);

} // namespace gstlrn