/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
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

  void setDmax(const VectorDouble &dmax) { _dmax = dmax; }
  void setFlagFill(bool flagFill) { _flagFill = flagFill; }
  void setFlagInter(bool flagInter) { _flagInter = flagInter; }
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

  virtual int _getNVar() const override;

private:
  int    _iattOut;
  VectorInt _iuids;
  int    _distType;
  VectorDouble _dmax;
  bool _flagFill;
  bool _flagInter;
  bool _flagLocate;
  ELoc _locatorType;
};

GSTLEARN_EXPORT int migrate(Db *dbin,
                            Db *dbout,
                            const String &name,
                            int dist_type = 1,
                            const VectorDouble &dmax = VectorDouble(),
                            bool flag_fill = false,
                            bool flag_inter = false,
                            const NamingConvention &namconv = NamingConvention(
                                "Migrate", false));
GSTLEARN_EXPORT int migrateMulti(Db *dbin,
                                 Db *dbout,
                                 const VectorString &names,
                                 int dist_type = 1,
                                 const VectorDouble &dmax = VectorDouble(),
                                 bool flag_fill = false,
                                 bool flag_inter = false,
                                 const NamingConvention &namconv = NamingConvention(
                                     "Migrate", false));
GSTLEARN_EXPORT int migrateByAttribute(Db *dbin,
                                       Db *dbout,
                                       const VectorInt &iatts = VectorInt(),
                                       int dist_type = 1,
                                       const VectorDouble &dmax = VectorDouble(),
                                       bool flag_fill = false,
                                       bool flag_inter = false,
                                       const NamingConvention &namconv = NamingConvention(
                                           "Migrate", false));
GSTLEARN_EXPORT int migrateByLocator(Db *dbin,
                                     Db *dbout,
                                     const ELoc &locatorType,
                                     int dist_type = 1,
                                     const VectorDouble &dmax = VectorDouble(),
                                     bool flag_fill = false,
                                     bool flag_inter = false,
                                     const NamingConvention &namconv = NamingConvention(
                                         "Migrate", false));
