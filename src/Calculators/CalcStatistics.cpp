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
#include "geoslib_old_f.h"

#include "Enum/ELoc.hpp"

#include "Basic/NamingConvention.hpp"
#include "Calculators/CalcStatistics.hpp"
#include "Calculators/ACalcDbToDb.hpp"
#include "Stats/Classical.hpp"
#include "Stats/Regression.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"

#include <math.h>

CalcStatistics::CalcStatistics()
    : ACalcDbToDb(),
      _iattOut(-1),
      _dboutMustBeGrid(false),
      _flagStats(false),
      _oper(EStatOption::UNKNOWN),
      _radius(0),
      _flagRegr(false),
      _flagCst(false),
      _regrMode(0),
      _nameResp(),
      _nameAux(),
      _model(nullptr)
{
}

CalcStatistics::~CalcStatistics()
{
}

bool CalcStatistics::_check()
{
  if (! ACalcDbToDb::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;

  int nvar = getDbin()->getLocatorNumber(ELoc::Z);
  if (nvar <= 0)
  {
    messerr("These methods require some variable to be defined");
    return false;
  }

  if (getDboutMustBeGrid())
  {
    if (! getDbout()->isGrid())
    {
      messerr("This method requires 'dbout' to be a Grid");
      return false;
    }
  }

  if (_flagRegr)
  {
    if (! _flagCst && _nameAux.empty())
    {
      messerr("This method requires Explanatory variables and/or constant term");
      return false;
    }
  }

  return true;
}

bool CalcStatistics::_preprocess()
{
  if (!ACalcDbToDb::_preprocess()) return false;
  
  if (_flagStats)
    _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, _getNVar(), 0.);

  if (_flagRegr)
    _iattOut = _addVariableDb(1, 1, ELoc::UNKNOWN, 0, 1, 0.);

  if (_iattOut < 0) return false;
  return true;
}

bool CalcStatistics::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  if (_flagStats)
    _renameVariable(2, VectorString(), ELoc::Z, _getNVar(), _iattOut, String(), 1);

  if (_flagRegr)
    _renameVariable(1, VectorString(), ELoc::Z, 1, _iattOut, String(), 1);
  return true;
}

void CalcStatistics::_rollback()
{
  _cleanVariableDb(1);
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcStatistics::_run()
{
  if (_flagStats)
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
    VectorString names = getDbin()->getNamesByLocator(ELoc::Z);
    if (dbStatisticsInGridTool(getDbin(), dbgrid, names, _oper, _radius, _iattOut))
      return false;
  }
  if (_flagRegr)
  {
    Regression reg = regression(getDbin(), _nameResp, _nameAux, _regrMode, _flagCst,
                                getDbout(), _model);
    if (reg.apply(getDbin(), _iattOut, _nameResp, _nameAux, _regrMode, _flagCst,
                  getDbout(), _model)) return false;
  }
  return true;
}

/****************************************************************************/
/*!
 **  Calculates the statistics on variables of an input Db per cell of an output Grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Input Db
 ** \param[in]  dbgrid  Output DbGrid
 ** \param[in]  oper    The statistical calculation
 ** \param[in]  radius  Neighborhood radius
 ** \param[in]  namconv Naming convention
 **
 *****************************************************************************/
int dbStatisticsOnGrid(Db *db,
                       DbGrid *dbgrid,
                       const EStatOption &oper,
                       int radius,
                       const NamingConvention &namconv)
{
  CalcStatistics stats;
  stats.setDbin(db);
  stats.setDbout(dbgrid);
  stats.setNamingConvention(namconv);

  stats.setFlagStats(true);
  stats.setDboutMustBeGrid(true);
  stats.setOper(oper);
  stats.setRadius(radius);

  // Run the calculator
  int error = (stats.run()) ? 0 : 1;
  return error;
}

int dbRegression(Db *db1,
                 const String& nameResp,
                 const VectorString& nameAux,
                 int mode,
                 bool flagCst,
                 Db *db2,
                 const Model* model,
                 const NamingConvention &namconv)
{
  if (db2 == nullptr) db2 = db1;

  CalcStatistics stats;
  stats.setDbin(db1);
  stats.setDbout(db2);
  stats.setNamingConvention(namconv);

  stats.setFlagRegr(true);
  stats.setRegrMode(mode);
  stats.setFlagCst(flagCst);
  stats.setName0(nameResp);
  stats.setNamaux(nameAux);
  stats.setModel(model);

  // Run the calculator
  int error = (stats.run()) ? 0 : 1;
  return error;
}
