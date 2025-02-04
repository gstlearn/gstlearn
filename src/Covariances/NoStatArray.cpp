
#include "Covariances/NoStatArray.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/VectorHelper.hpp"
#include "geoslib_define.h"

NoStatArray::NoStatArray(const Db *dbref,const String& colname)
:_dbNoStat(dbref)
,_colName(colname)
{}


String NoStatArray::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << ANoStat::toString(strfmt);
  AStringFormat sf;
  if (strfmt != nullptr) sf = *strfmt;
  return sstr.str();
}


void NoStatArray::_informField(const VectorVectorDouble& coords,
                              VectorDouble& tab, bool verbose)
{
  // Identify the attribute in the Db

  int iatt = _dbNoStat->getUID(_colName);
  if (iatt < 0)
  {
    messerr("The Non-stationary attribute  is not defined in _dbNoStat anymore");
    return;
  }

  // Migrate the information from Db onto the Vertex locations

  if (_dbNoStat->isGrid())
  {
    const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(_dbNoStat);
    if (migrateGridToCoor(dbgrid, iatt, coords, tab)) return;
  }
  else
  {
    if (expandPointToCoor(_dbNoStat, iatt, coords, tab)) return;
  }

  int ndef = VH::countUndefined(tab);
  if (ndef > 0)
  {

    // Calculate local statistics

    double mean = VH::mean(tab);
    if (FFFF(mean))
    {
      messerr("This Non-Stationary parameter is not valid");
      return;
    }

    if (verbose)
    {
      message("For Non-Stationary Parameter, there are %d undefined values\n",
              ndef);
      message("They have been replaced by its average value (%lf)\n", mean);
    }

    // Modify the TEST values to the mean value

    VH::fillUndef(tab, mean);
  }

  // Printout some statistics (optional)

  if (verbose)
  {
    char str[LONG_SIZE];
    (void) gslSPrintf(str,
                      "Statistics for Non-Stationary Parameter on Mesh");
    VH::dumpStats(str,tab);
  }

}
