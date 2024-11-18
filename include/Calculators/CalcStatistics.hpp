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
class Model;

class GSTLEARN_EXPORT CalcStatistics: public ACalcDbToDb
{
public:
  CalcStatistics();
  CalcStatistics(const CalcStatistics &r) = delete;
  CalcStatistics& operator=(const CalcStatistics &r) = delete;
  virtual ~CalcStatistics();

  bool getDboutMustBeGrid() const { return _dboutMustBeGrid; }
  void setDboutMustBeGrid(bool dboutMustBeGrid) { _dboutMustBeGrid = dboutMustBeGrid; }

  void setFlagStats(bool flagStats) { _flagStats = flagStats; }
  void setRadius(int radius) { _radius = radius; }
  void setOper(const EStatOption &oper) { _oper = oper; }

  void setFlagRegr(bool flagRegr) { _flagRegr = flagRegr; }
  void setFlagCst(bool flagCst) { _flagCst = flagCst; }
  void setName0(const String &name0) { _nameResp = name0; }
  void setNamaux(const VectorString &namaux) { _nameAux = namaux; }
  void setRegrMode(int regrMode) { _regrMode = regrMode; }
  void setModel(const Model* model) { _model = model; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

private:
  int    _iattOut;
  bool   _dboutMustBeGrid;
  bool   _flagStats;
  EStatOption _oper;
  int    _radius;
  bool   _flagRegr;
  bool   _flagCst;
  int    _regrMode;
  String _nameResp;
  VectorString _nameAux;
  const Model* _model;
};

GSTLEARN_EXPORT int dbStatisticsOnGrid(Db *db,
                                       DbGrid *dbgrid,
                                       const EStatOption &oper,
                                       int radius = 0,
                                       const NamingConvention &namconv = NamingConvention("Stats"));
GSTLEARN_EXPORT int dbRegression(Db *db1,
                                 const String& nameResp,
                                 const VectorString& nameAux,
                                 int mode = 0,
                                 bool flagCst = true,
                                 Db *db2 = nullptr,
                                 const Model* model = nullptr,
                                 const NamingConvention &namconv = NamingConvention("Regr"));
