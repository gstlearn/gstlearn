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
class Model;
class GSTLEARN_EXPORT CalcStatistics: public ACalcDbToDb
{
public:
  CalcStatistics();
  CalcStatistics(const CalcStatistics& r)            = delete;
  CalcStatistics& operator=(const CalcStatistics& r) = delete;
  virtual ~CalcStatistics();

  bool getDboutMustBeGrid() const { return _dboutMustBeGrid; }
  void setDboutMustBeGrid(bool dboutMustBeGrid) { _dboutMustBeGrid = dboutMustBeGrid; }

  void setFlagStats(bool flagStats) { _flagStats = flagStats; }
  void setRadius(Id radius) { _radius = radius; }
  void setOper(const EStatOption& oper) { _oper = oper; }

  void setFlagRegr(bool flagRegr) { _flagRegr = flagRegr; }
  void setFlagCst(bool flagCst) { _flagCst = flagCst; }
  void setName0(const String& name0) { _nameResp = name0; }
  void setNamaux(const VectorString& namaux) { _nameAux = namaux; }
  void setRegrMode(Id regrMode) { _regrMode = regrMode; }
  void setModel(const Model* model) { _model = model; }

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

private:
  Id _iattOut;
  bool _dboutMustBeGrid;
  bool _flagStats;
  EStatOption _oper;
  Id _radius;
  bool _flagRegr;
  bool _flagCst;
  Id _regrMode;
  String _nameResp;
  VectorString _nameAux;
  const Model* _model;
};

GSTLEARN_EXPORT Id dbStatisticsOnGrid(Db* db,
                                       DbGrid* dbgrid,
                                       const EStatOption& oper,
                                       Id radius                      = 0,
                                       const NamingConvention& namconv = NamingConvention("Stats"));
GSTLEARN_EXPORT Id dbRegression(Db* db1,
                                 const String& nameResp,
                                 const VectorString& nameAux,
                                 Id mode                        = 0,
                                 bool flagCst                    = true,
                                 Db* db2                         = nullptr,
                                 const Model* model              = nullptr,
                                 const NamingConvention& namconv = NamingConvention("Regr"));
} // namespace gstlrn