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
#include "OutputFormat/GridEclipse.hpp"
#include "OutputFormat/AOF.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>

GridEclipse::GridEclipse(const char* filename, const Db* db)
  : AOF(filename, db)
{
}

GridEclipse::GridEclipse(const GridEclipse& r)
    : AOF(r)
{
}

GridEclipse& GridEclipse::operator=(const GridEclipse& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
  }
  return *this;
}

GridEclipse::~GridEclipse()
{
}

int GridEclipse::writeInFile()
{
  static int nbyline = 6;
  static double valtest = -9999.;

  /* Open the file */

  if (_fileWriteOpen()) return 1;

  // Preliminary calculations

  int nxyz = 1;
  for (int idim = 0; idim < _dbgrid->getNDim(); idim++)
    nxyz *= _dbgrid->getNX(idim);

  /* Write a comment */

  fprintf(_file, "Facies\n");

  /* Write the set of values */

  int ninline = 0;
  for (int i = 0; i < nxyz; i++)
  {
    double valprt = valtest;
    if (_dbgrid->getSelection(i))
    {
      double value = _dbgrid->getArray(i, _cols[0]);
      if (!FFFF(value)) valprt = value;
    }
    fprintf(_file, "%lf ", valprt);
    ninline++;
    if (ninline == nbyline)
    {
      fprintf(_file, "\n");
      ninline = 0;
    }
  }
  if (ninline > 0) fprintf(_file, "\n");

  _fileClose();
  return 0;
}
