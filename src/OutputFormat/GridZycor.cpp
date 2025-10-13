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
#include "OutputFormat/GridZycor.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "Core/io.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Enum/ELoadBy.hpp"
#include "OutputFormat/AOF.hpp"

#include <cstdio>
#include <cstring>

#define ZYCOR_NULL_CH "  0.1000000E+31"

namespace gstlrn
{
GridZycor::GridZycor(const char* filename, const Db* db)
  : AOF(filename, db)
{
}

GridZycor::GridZycor(const GridZycor& r)
  : AOF(r)
{
}

GridZycor& GridZycor::operator=(const GridZycor& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
  }
  return *this;
}

GridZycor::~GridZycor()
{
}

Id GridZycor::writeInFile()
{
  Id nx[2];
  double rbid, x0[2], xf[2], dx[2];
  double buff[5]; /* Size = nbyline */
  Id card_max = 100;
  char card[100]; /* Size = nbyline * 20 */
  static Id nbyline     = 5;
  static double testval = MAXIMUM_BIG;

  /* Open the file */

  if (_fileWriteOpen()) return 1;

  /* Write a comment */

  fprintf(_file, "!\n");
  fprintf(_file, "!  File created by gstlearn package\n");
  fprintf(_file, "!\n");

  /* Title line */

  fprintf(_file, "@GRID ZYCOR FILE    ,   GRID,  %ld\n", nbyline);
  fprintf(_file, "     15, %13lg,    ,    0,     1\n", testval);

  /* Grid description */

  for (Id i = 0; i < 2; i++)
  {
    nx[i] = _dbgrid->getNX(i);
    x0[i] = _dbgrid->getX0(i);
    dx[i] = _dbgrid->getDX(i);
    xf[i] = x0[i] + (nx[i] - 1) * dx[i];
  }

  rbid = 0.;
  fprintf(_file, "%6ld, %6ld, %13lf, %13lf, %13lf, %13lf\n", nx[1], nx[0], x0[0],
          xf[0], x0[1], xf[1]);
  fprintf(_file, " %15lf, %15lf, %15lf\n", rbid, rbid, rbid);
  fprintf(_file, "@\n");

  /* The set of values */

  for (Id jj = nx[0] - 1; jj >= 0; jj--)
  {
    Id kk = 0;
    Id ii = ((nx[1] * nx[0]) - (jj + 1));
    for (Id loop = 1; loop <= nx[1]; loop++)
    {
      buff[kk++] = _dbgrid->getArray(ii, _cols[0]);
      ii -= nx[0];
      if (kk == nbyline)
      {
        for (Id yy = 0; yy < nbyline; yy++)
        {
          Id ind = yy * 15;
          if (!FFFF(buff[yy]))
          {
            gslSPrintf(&card[ind], card_max, "%15g", buff[yy]);
          }
          else
          {
            gslStrcpy(&card[ind], card_max, ZYCOR_NULL_CH);
          }
        }
        gslSPrintf(&card[15 * nbyline], card_max, "\n");
        fprintf(_file, "%s", card);
        kk = 0;
      }
    }

    if (kk > 0)
    {
      for (Id yy = 0; yy < kk; yy++)
      {
        Id ind = yy * 15;
        if (!FFFF(buff[yy]))
        {
          gslSPrintf(&card[ind], card_max, "%15g", buff[yy]);
        }
        else
        {
          gslStrcpy(&card[ind], card_max, ZYCOR_NULL_CH);
        }
      }
      gslSPrintf(&card[15 * kk], card_max, "\n");
      fprintf(_file, "%s", card);
    }
  }

  _fileClose();
  return 0;
}

DbGrid* GridZycor::readGridFromFile()
{
  DbGrid* dbgrid = nullptr;

  /* Open the file */

  if (_fileReadOpen()) return dbgrid;

  /* Define the delimitors */

  _token_delimitors('!', ',', ' ');

  /* Read the lines */

  String string;
  if (_record_read(_file, "%s", &string)) return dbgrid;
  if (string[0] != '@')
  {
    messerr("Missing string starting with (@). Instead: '%s'", string.data());
    return dbgrid;
  }
  if (_record_read(_file, "%s", &string)) return dbgrid;
  if (string != "GRID")
  {
    messerr("Missing string (GRID). Instead: '%s'", string.data());
    return dbgrid;
  }

  double rbid1, rbid2, rbid3, test, value;
  Id nval, ibid1, ibid2, ibid3;
  VectorInt nx(2);
  VectorDouble x0(2);
  VectorDouble xf(2);
  if (_record_read(_file, "%ld", &nval)) return dbgrid;
  if (_record_read(_file, "%ld", &ibid1)) return dbgrid;
  if (_record_read(_file, "%lg", &test)) return dbgrid;
  if (_record_read(_file, "%s", &string)) return dbgrid;
  if (_record_read(_file, "%ld", &ibid2)) return dbgrid;
  if (_record_read(_file, "%ld", &ibid3)) return dbgrid;
  if (_record_read(_file, "%ld", nx.data() + 1)) return dbgrid;
  if (_record_read(_file, "%ld", nx.data() + 0)) return dbgrid;
  if (_record_read(_file, "%lf", x0.data() + 0)) return dbgrid;
  if (_record_read(_file, "%lf", xf.data() + 0)) return dbgrid;
  if (_record_read(_file, "%lf", x0.data() + 1)) return dbgrid;
  if (_record_read(_file, "%lf", xf.data() + 1)) return dbgrid;
  if (_record_read(_file, "%lf", &rbid1)) return dbgrid;
  if (_record_read(_file, "%lf", &rbid2)) return dbgrid;
  if (_record_read(_file, "%lf", &rbid3)) return dbgrid;

  if (_record_read(_file, "%s", &string)) return dbgrid;
  if (string != "@")
  {
    messerr("Missing string (@). Instead: %s", string.data());
    return dbgrid;
  }

  /* Final calculations */

  VectorDouble dx(2);
  dx[0] = (xf[0] - x0[0]) / (nx[0] - 1);
  dx[1] = (xf[1] - x0[1]) / (nx[1] - 1);

  /* Reset the delimitors */

  _token_delimitors('#', ' ', ' ');

  /* Core allocation */

  Id size = nx[0] * nx[1];
  VectorDouble tab(size);

  /* Read the array of real values */

  for (Id ix = 0; ix < nx[0]; ix++)
    for (Id iy = 0; iy < nx[1]; iy++)
    {
      if (_record_read(_file, "%lf", &value)) break;
      if (value == test) value = TEST;
      tab[(nx[1] - iy - 1) * nx[0] + ix] = value;
    }

  dbgrid = new DbGrid();
  dbgrid->reset(nx, dx, x0, VectorDouble(), ELoadBy::SAMPLE, tab);

  // Close the file

  _fileClose();

  return dbgrid;
}
} // namespace gstlrn
