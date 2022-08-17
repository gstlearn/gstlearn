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
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "Basic/NamingConvention.hpp"
#include "Calculators/ACalcDb2Db.hpp"
#include "Calculators/CalcStatistics.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Db/ELoc.hpp"

#include <math.h>

CalcStatistics::CalcStatistics()
    : ACalcDb2Db(),
      _iattOut(-1)
{
}

CalcStatistics::~CalcStatistics()
{
}

bool CalcStatistics::_check()
{
  if (! ACalcDb2Db::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;

  if (! getDbout()->isGrid())
  {
    messerr("This method requires 'dbout' to be a Grid");
    return false;
  }
  if (_getNVar() <= 0)
  {
    messerr("These methods require some variable to be defined");
    return false;
  }
  return true;
}

bool CalcStatistics::_preprocess()
{
  int nvar = _getNVar();
  _iattOut = _addVariableDb(2, 1, ELoc::UNKNOWN, nvar, 0.);
  if (_iattOut < 0) return false;
  return true;
}

bool CalcStatistics::_postprocess()
{
  _renameVariable(1, _iattOut, String(), 1);
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
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  calc_stats_grid(getDbin(), dbgrid, _oper, _iattOut, _radius);

  return true;
}

int dbStatisticsOnGrid(Db *db,
                        DbGrid *dbgrid,
                        const EStatOption &oper,
                        int radius,
                        const NamingConvention& namconv)
{
  CalcStatistics stats;
  stats.setDbin(db);
  stats.setDbout(dbgrid);
  stats.setNamingConvention(namconv);

  stats.setOper(oper);
  stats.setRadius(radius);

  // Run the calculator
  int error = (stats.run()) ? 0 : 1;
  return error;
}
