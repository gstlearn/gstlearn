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


class Db;
class DbGrid;
class EStatOption;

class GSTLEARN_EXPORT CalcMigrate: public ACalcDbToDb
{
public:
  CalcMigrate();
  CalcMigrate(const CalcMigrate &r) = delete;
  CalcMigrate& operator=(const CalcMigrate &r) = delete;
  virtual ~CalcMigrate();

  void setDmax(const VectorDouble& dmax) { _dmax = dmax; }
  void setFlagFill(bool flagFill) { _flagFill = flagFill; }
  void setFlagInter(bool flagInter) { _flagInter = flagInter; }
  void setFlagBall(bool flagBall) { _flagBall = flagBall; }
  void setDistType(int dist_type) { _distType = dist_type; }
  void setIuids(const VectorInt &iuids) { _iuids = iuids; }
  void setFlagLocate(bool flagLocate) { _flagLocate = flagLocate; }
  void setLocatorType(const ELoc &locatorType) { _locatorType = locatorType; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  static int _migrate(Db* db1,
                      Db* db2,
                      int iatt1,
                      int iatt2,
                      int distType,
                      const VectorDouble& dmax,
                      bool flag_fill,
                      bool flag_inter,
                      bool flag_ball);
  static int _migratePointToGrid(Db* db_point,
                                 DbGrid* db_grid,
                                 int iatt,
                                 int distType,
                                 const VectorDouble& dmax,
                                 VectorDouble& tab);
  static int _expandPointToPointBall(Db* db1,
                                     Db* db2,
                                     int iatt,
                                     int distType,
                                     const VectorDouble& dmax,
                                     VectorDouble& tab);
  static int _migrateGridToGrid(DbGrid* db_gridin,
                                DbGrid* db_gridout,
                                int iatt,
                                int distType,
                                const VectorDouble& dmax,
                                VectorDouble& tab);
  static int _expandPointToPoint(Db* db1,
                                 Db* db2,
                                 int iatt,
                                 int distType,
                                 const VectorDouble& dmax,
                                 VectorDouble& tab);
  static int _expandGridToGrid(DbGrid* db_gridin,
                               DbGrid* db_gridout,
                               int iatt,
                               int distType,
                               const VectorDouble& dmax,
                               VectorDouble& tab);
  static int _interpolateGridToPoint(DbGrid* db_grid,
                                     Db* db_point,
                                     int iatt,
                                     int distType,
                                     const VectorDouble& dmax,
                                     VectorDouble& tab);
  static int _migrateGridToPoint(DbGrid* db_grid,
                                 Db* db_point,
                                 int iatt,
                                 int distType,
                                 const VectorDouble& dmax,
                                 VectorDouble& tab);

private:
  int    _iattOut;
  VectorInt _iuids;
  int    _distType;
  VectorDouble _dmax;
  bool _flagFill;
  bool _flagInter;
  bool _flagLocate;
  bool _flagBall;
  ELoc _locatorType;
};

GSTLEARN_EXPORT int migrate(Db *dbin,
                            Db *dbout,
                            const String &name,
                            int dist_type = 1,
                            const VectorDouble& dmax = VectorDouble(),
                            bool flag_fill = false,
                            bool flag_inter = false,
                            bool flag_ball = false,
                            const NamingConvention &namconv = NamingConvention(
                                "Migrate", false));
GSTLEARN_EXPORT int migrateMulti(Db *dbin,
                                 Db *dbout,
                                 const VectorString &names,
                                 int dist_type = 1,
                                 const VectorDouble& dmax = VectorDouble(),
                                 bool flag_fill = false,
                                 bool flag_inter = false,
                                 bool flag_ball = false,
                                 const NamingConvention &namconv = NamingConvention(
                                     "Migrate"));
GSTLEARN_EXPORT int migrateByAttribute(Db *dbin,
                                       Db *dbout,
                                       const VectorInt &iatts = VectorInt(),
                                       int dist_type = 1,
                                       const VectorDouble& dmax = VectorDouble(),
                                       bool flag_fill = false,
                                       bool flag_inter = false,
                                       bool flag_ball = false,
                                       const NamingConvention &namconv = NamingConvention(
                                           "Migrate"));
GSTLEARN_EXPORT int migrateByLocator(Db *dbin,
                                     Db *dbout,
                                     const ELoc &locatorType,
                                     int dist_type = 1,
                                     const VectorDouble& dmax = VectorDouble(),
                                     bool flag_fill = false,
                                     bool flag_inter = false,
                                     bool flag_ball = false,
                                     const NamingConvention &namconv = NamingConvention(
                                         "Migrate"));
GSTLEARN_EXPORT int manageExternalInformation(int mode,
                                              const ELoc &locatorType,
                                              Db *dbin,
                                              Db *dbout,
                                              bool* flag_created);
GSTLEARN_EXPORT int interpolateVariableToPoint(DbGrid *db_grid,
                                               int iatt,
                                               int np,
                                               const double *xp,
                                               const double *yp,
                                               const double *zp,
                                               double *tab);
GSTLEARN_EXPORT double* dbgridLineSampling(DbGrid *dbgrid,
                                       const double *x1,
                                       const double *x2,
                                       int ndisc,
                                       int ncut,
                                       const double *cuts,
                                       int *nval_ret);
GSTLEARN_EXPORT int expandPointToGrid(Db *db_point,
                                      DbGrid *db_grid,
                                      int iatt,
                                      int iatt_time,
                                      int iatt_angle,
                                      int iatt_scaleu,
                                      int iatt_scalev,
                                      int iatt_scalew,
                                      int flag_index,
                                      int distType,
                                      const VectorDouble &dmax,
                                      VectorDouble &tab);
GSTLEARN_EXPORT int pointToBlock(Db *dbpoint,
                                 DbGrid *dbgrid,
                                 int option,
                                 int flag_size,
                                 int iatt_time,
                                 int iatt_size,
                                 int iatt_angle,
                                 int iatt_scaleu,
                                 int iatt_scalev,
                                 int iatt_scalew);
GSTLEARN_EXPORT int migrateGridToCoor(const DbGrid *db_grid,
                                      int iatt,
                                      const VectorVectorDouble &coords,
                                      VectorDouble &tab);
GSTLEARN_EXPORT int expandPointToCoor(const Db *db1,
                                      int iatt,
                                      const VectorVectorDouble &coords,
                                      VectorDouble &tab);

