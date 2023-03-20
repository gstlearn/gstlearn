/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "OutputFormat/GridXYZ.hpp"
#include "OutputFormat/AOF.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>

GridXYZ::GridXYZ(const char* filename, const Db* db)
  : AOF(filename, db)
{
}

GridXYZ::GridXYZ(const GridXYZ& r)
    : AOF(r)
{
}

GridXYZ& GridXYZ::operator=(const GridXYZ& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
  }
  return *this;
}

GridXYZ::~GridXYZ()
{
}

int GridXYZ::writeInFile()
{
  /* Open the file */

  if (_fileWriteOpen()) return 1;

  /* Write a comment */

  fprintf(_file, "FDASCII 0 0 0 0 1E30\n");
  fprintf(_file, "->\n");

  /* Write the set of values */

  int lec = 0;
  for (int ix = 0; ix < _dbgrid->getNX(0); ix++)
    for (int iy = 0; iy < _dbgrid->getNX(1); iy++)
    {
      for (int i = 0; i < _dbgrid->getNDim(); i++)
        fprintf(_file, "%lf,", _dbgrid->getCoordinate(lec, i));
      double value = _dbgrid->getArray(lec, _cols[0]);
      if (FFFF(value))
        fprintf(_file, "1E+30\n");
      else
        fprintf(_file, "%lf\n", value);
      lec++;
    }

  // Close the file

  _fileClose();
  return 0;
}
