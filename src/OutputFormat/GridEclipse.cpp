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
#include "OutputFormat/GridEclipse.hpp"
#include "OutputFormat/AOF.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

namespace gstlrn
{
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

Id GridEclipse::writeInFile()
{
  static Id nbyline = 6;
  static double valtest = -9999.;

  /* Open the file */

  if (_fileWriteOpen()) return 1;

  // Preliminary calculations

  Id nxyz = 1;
  for (Id idim = 0; idim < _dbgrid->getNDim(); idim++)
    nxyz *= _dbgrid->getNX(idim);

  /* Write a comment */

  fprintf(_file, "Facies\n");

  /* Write the set of values */

  Id ninline = 0;
  for (Id i = 0; i < nxyz; i++)
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
}