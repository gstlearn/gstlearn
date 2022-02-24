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
#include "OutputFormat/GridZycor.hpp"
#include "OutputFormat/AOF.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>

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

int GridZycor::dumpFile()
{
  int nx[2];
  double rbid, x0[2], xf[2], dx[2];
  double buff[5]; /* Size = nbyline */
  char card[100]; /* Size = nbyline * 20 */
  static int nbyline = 5;
  static double testval = 1.e30;

  /* Open the file */

  if (_fileOpen()) return 1;

  /* Write a comment */

  fprintf(_file, "!\n");
  fprintf(_file, "!  File created by RGeostats package\n");
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
