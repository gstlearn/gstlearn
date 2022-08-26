/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
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
  void setLdmax(int ldmax) { _ldmax = ldmax; }
  void setIuids(const VectorInt &iuids) { _iuids = iuids; }

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
  int    _ldmax;
  VectorDouble _dmax;
  bool _flagFill;
  bool _flagInter;
};

GSTLEARN_EXPORT int migrate(Db *db1,
                            Db *db2,
                            const String &name,
                            int ldmax = 1,
                            const VectorDouble &dmax = VectorDouble(),
                            int flag_fill = 0,
                            int flag_inter = 0,
                            const NamingConvention& namconv = NamingConvention("Migrate"));
GSTLEARN_EXPORT int migrateByAttribute(Db *db1,
                                       Db *db2,
                                       const VectorInt &iatts = VectorInt(),
                                       int ldmax = 1,
                                       const VectorDouble &dmax = VectorDouble(),
                                       int flag_fill = false,
                                       int flag_inter = false,
                                       const NamingConvention& namconv = NamingConvention("Migrate"));
GSTLEARN_EXPORT int migrateByLocator(Db *db1,
                                     Db *db2,
                                     const ELoc &locatorType,
                                     int ldmax = 1,
                                     const VectorDouble &dmax = VectorDouble(),
                                     int flag_fill = false,
                                     int flag_inter = false,
                                     const NamingConvention& namconv = NamingConvention("Migrate"));
