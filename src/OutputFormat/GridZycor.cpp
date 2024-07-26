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
#include "geoslib_f_private.h"

#include "Enum/ELoadBy.hpp"

#include "OutputFormat/GridZycor.hpp"
#include "OutputFormat/AOF.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>
#include <stdio.h>

#define ZYCOR_NULL_CH "  0.1000000E+31"

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

int GridZycor::writeInFile()
{
  int nx[2];
  double rbid, x0[2], xf[2], dx[2];
  double buff[5]; /* Size = nbyline */
  char card[100]; /* Size = nbyline * 20 */
  static int nbyline = 5;
  static double testval = 1.e30;

  /* Open the file */

  if (_fileWriteOpen()) return 1;

  /* Write a comment */

  fprintf(_file, "!\n");
  fprintf(_file, "!  File created by gstlearn package\n");
  fprintf(_file, "!\n");

  /* Title line */

  fprintf(_file, "@GRID ZYCOR FILE    ,   GRID,  %d\n", nbyline);
  fprintf(_file, "     15, %13lg,    ,    0,     1\n", testval);

  /* Grid description */

  for (int i = 0; i < 2; i++)
  {
    nx[i] = _dbgrid->getNX(i);
    x0[i] = _dbgrid->getX0(i);
    dx[i] = _dbgrid->getDX(i);
    xf[i] = x0[i] + (nx[i] - 1) * dx[i];
  }

  rbid = 0.;
  fprintf(_file, "%6d, %6d, %13lf, %13lf, %13lf, %13lf\n", nx[1], nx[0], x0[0],
          xf[0], x0[1], xf[1]);
  fprintf(_file, " %15lf, %15lf, %15lf\n", rbid, rbid, rbid);
  fprintf(_file, "@\n");

  /* The set of values */

  for (int jj = nx[0] - 1; jj >= 0; jj--)
  {
    int kk = 0;
    int ii = ((nx[1] * nx[0]) - (jj + 1));
    for (int loop = 1; loop <= nx[1]; loop++)
    {
      buff[kk++] = _dbgrid->getArray(ii, _cols[0]);
      ii -= nx[0];
      if (kk == nbyline)
      {
        for (int yy = 0; yy < nbyline; yy++)
        {
          int ind = yy * 15;
          if (!FFFF(buff[yy]))
          {
            gslSPrintf(&card[ind], "%15g", buff[yy]);
          }
          else
          {
            memcpy(&card[ind], (char*) ZYCOR_NULL_CH, 15);
          }
        }
        gslSPrintf(&card[15 * nbyline], "\n");
        fprintf(_file, "%s", card);
        kk = 0;
      }
    }

    if (kk > 0)
    {
      for (int yy = 0; yy < kk; yy++)
      {
        int ind = yy * 15;
        if (!FFFF(buff[yy]))
        {
          gslSPrintf(&card[ind], "%15g", buff[yy]);
        }
        else
        {
          memcpy(&card[ind], (char*) ZYCOR_NULL_CH, 15);
        }
      }
      gslSPrintf(&card[15 * kk], "\n");
      fprintf(_file, "%s", card);
    }
  }

  _fileClose();
  return 0;
}

DbGrid*  GridZycor::readGridFromFile()
{
  DbGrid* dbgrid = nullptr;
  char string[100];
  double xf[2], rbid1, rbid2, rbid3, test, value;
  int nval, ibid1, ibid2, ibid3;
  VectorInt nx(2);
  VectorDouble dx(2);
  VectorDouble x0(2);

  /* Open the file */

  if (_fileReadOpen()) return dbgrid;

   /* Define the delimitors */

   _file_delimitors('!', ',', '_');

   /* Read the lines */

   if (_record_read(_file, "%s", string)) return dbgrid;
   if (string[0] != '@')
   {
     messerr("Missing string starting with (@). Instead: '%s'", string);
     return dbgrid;
   }
   if (_record_read(_file, "%s", string)) return dbgrid;
   if (strcmp(string, "GRID") != 0)
   {
     messerr("Missing string (GRID). Instead: '%s'", string);
     return dbgrid;
   }
   if (_record_read(_file, "%d", &nval)) return dbgrid;
   if (_record_read(_file, "%d", &ibid1)) return dbgrid;
   if (_record_read(_file, "%lg", &test)) return dbgrid;
   if (_record_read(_file, "%s", string)) return dbgrid;
   if (_record_read(_file, "%d", &ibid2)) return dbgrid;
   if (_record_read(_file, "%d", &ibid3)) return dbgrid;
   if (_record_read(_file, "%d", &nx[1])) return dbgrid;
   if (_record_read(_file, "%d", &nx[0])) return dbgrid;
   if (_record_read(_file, "%lf", &x0[0])) return dbgrid;
   if (_record_read(_file, "%lf", &xf[0])) return dbgrid;
   if (_record_read(_file, "%lf", &x0[1])) return dbgrid;
   if (_record_read(_file, "%lf", &xf[1])) return dbgrid;
   if (_record_read(_file, "%lf", &rbid1)) return dbgrid;
   if (_record_read(_file, "%lf", &rbid2)) return dbgrid;
   if (_record_read(_file, "%lf", &rbid3)) return dbgrid;

   if (_record_read(_file, "%s", string)) return dbgrid;
   if (strcmp(string, "@") != 0)
   {
     messerr("Missing string (@). Instead: %s", string);
     return dbgrid;
   }

   /* Final calculations */

   dx[0] = (xf[0] - x0[0]) / (nx[0] - 1);
   dx[1] = (xf[1] - x0[1]) / (nx[1] - 1);

   /* Reset the delimitors */

   _file_delimitors('#', ' ', ' ');

   /* Core allocation */

   int size = nx[0] * nx[1];
   VectorDouble tab(size);

   /* Read the array of real values */

   for (int ix = 0; ix < nx[0]; ix++)
     for (int iy = 0; iy < nx[1]; iy++)
     {
       if (_record_read(_file, "%lf", &value)) break;
       if (value == test) value = TEST;
       tab[(nx[1] - iy - 1) * nx[0] + ix] = value;
     }

   dbgrid = new DbGrid();
   dbgrid->reset(nx,dx,x0,VectorDouble(),ELoadBy::SAMPLE,tab);

  // Close the file

  _fileClose();

  return dbgrid;
}
