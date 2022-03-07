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
#include "OutputFormat/GridArcGis.hpp"
#include "OutputFormat/AOF.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>

GridArcGis::GridArcGis(const char* filename, const Db* db)
  : AOF(filename, db)
{
}

GridArcGis::GridArcGis(const GridArcGis& r)
    : AOF(r)
{
}

GridArcGis& GridArcGis::operator=(const GridArcGis& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
  }
  return *this;
}

GridArcGis::~GridArcGis()
{
}

bool GridArcGis::isAuthorized() const
{
  if (_dbgrid->getDX(0) != _dbgrid->getDX(1))
  {
    messerr("This function requires the Grid Mesh to be square");
    return false;
  }
  return true;
}

int GridArcGis::writeInFile()
{
  static double noValue = -9999.;

  /* Open the file */

  if (_fileWriteOpen()) return 1;

  /* Write a comment */

  fprintf(_file, "NCOLS %d\n", _dbgrid->getNX(0));
  fprintf(_file, "NROWS %d\n", _dbgrid->getNX(1));
  fprintf(_file, "XLLCORNER %lf\n", _dbgrid->getX0(0));
  fprintf(_file, "YLLCORNER %lf\n", _dbgrid->getX0(1));
  fprintf(_file, "CELLSIZE %lf\n", _dbgrid->getDX(0));
  fprintf(_file, "NODATA_VALUE %lf\n", noValue);

  /* Write the set of values */

  int lec = 0;
  for (int ix = 0; ix < _dbgrid->getNX(0); ix++)
    for (int iy = 0; iy < _dbgrid->getNX(1); iy++)
    {
      double value = _dbgrid->getArray(lec, _cols[0]);
      if (FFFF(value))
        fprintf(_file, "%lf\n", noValue);
      else
        fprintf(_file, "%lf\n", value);
      lec++;
    }

  /* Close the file */

  _fileClose();
  return 0;
}

