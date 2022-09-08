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
#include "geoslib_old_f.h"

#include "Basic/NamingConvention.hpp"
#include "Calculators/CalcStatistics.hpp"
#include "Calculators/ACalcDbToDb.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Db/ELoc.hpp"

#include <math.h>

CalcStatistics::CalcStatistics()
    : ACalcDbToDb(),
      _iattOut(-1),
      _flagStats(false),
      _oper(EStatOption::UNKNOWN),
      _radius(0),
      _flagRegr(false),
      _flagCste(false),
      _regrMode(0),
      _namaux()
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

  if (_getNVar() <= 0)
  {
    messerr("These methods require some variable to be defined");
    return false;
  }

  if (_flagStats)
  {
    if (! getDbout()->isGrid())
    {
      messerr("This method requires 'dbout' to be a Grid");
      return false;
    }
  }
  if (_flagRegr)
  {
    if (_getNVar() != 1)
    {
      messerr("This method requires a single variable in 'Dbin'");
      return false;
    }
    if (! _flagCste && _namaux.empty())
    {
      messerr("This method requires Explanatory variables and/or constant term");
      return false;
    }
  }

  return true;
}

bool CalcStatistics::_preprocess()
{
  int nvar = _getNVar();

  if (_flagStats)
    _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, nvar, 0.);

  if (_flagRegr)
    _iattOut = _addVariableDb(1, 1, ELoc::UNKNOWN, 0, nvar, 0.);

  if (_iattOut < 0) return false;
  return true;
}

bool CalcStatistics::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  if (_flagStats)
    _renameVariable(2, 1, _iattOut, String(), 1);

  if (_flagRegr)
    _renameVariable(1, 1, _iattOut, String(), 1);
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
    calc_stats_grid(getDbin(), dbgrid, _oper, _radius, _iattOut);
  }

  if (_flagRegr)
  {
    int icol0 = getDbin()->getUIDByLocator(ELoc::Z, 0);
    VectorInt icols = getDbout()->getUIDs(_namaux);
    calc_regression(getDbin(), getDbout(), _regrMode, icol0, icols, _flagCste, _iattOut);
  }
  return true;
}

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
  stats.setOper(oper);
  stats.setRadius(radius);

  // Run the calculator
  int error = (stats.run()) ? 0 : 1;
  return error;
}

int dbRegression(Db *db1,
                 Db *db2,
                 int mode,
                 const VectorString &namaux,
                 bool flagCste,
                 const NamingConvention &namconv)
{
  CalcStatistics stats;
  stats.setDbin(db1);
  stats.setDbout(db2);
  stats.setNamingConvention(namconv);

  stats.setFlagRegr(true);
  stats.setRegrMode(mode);
  stats.setFlagCste(flagCste);
  stats.setNamaux(namaux);

  // Run the calculator
  int error = (stats.run()) ? 0 : 1;
  return error;
}

