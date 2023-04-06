/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "OutputFormat/GridIrap.hpp"
#include "OutputFormat/AOF.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>

#define N_SAMPLE(nx,nsample) ((int) ((nx-1) / nsample) + 1)

GridIrap::GridIrap(const char* filename, const Db* db)
  : AOF(filename, db)
  , _nsamplex(1)
  , _nsampley(1)
{
}

GridIrap::GridIrap(const GridIrap& r)
    : AOF(r),
      _nsamplex(r._nsamplex),
      _nsampley(r._nsampley)
{
}

GridIrap& GridIrap::operator=(const GridIrap& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
    _nsamplex = r._nsamplex;
    _nsampley = r._nsampley;
  }
  return *this;
}

GridIrap::~GridIrap()
{
}

int GridIrap::writeInFile()
{
  VectorInt indg(2);

  /* Open the file */

  if (_fileWriteOpen()) return 1;

  /* Preliminary calculations */

  int nx = N_SAMPLE(_dbgrid->getNX(0), _nsamplex);
  int ny = N_SAMPLE(_dbgrid->getNX(1), _nsampley);
  double dx = _dbgrid->getDX(0) * _nsamplex;
  double dy = _dbgrid->getDX(1) * _nsampley;
  double xmin = _dbgrid->getX0(0);
  double ymin = _dbgrid->getX0(1);
  double xmax = xmin + dx * (nx - 1);
  double ymax = ymin + dy * (ny - 1);

  /* Write the header */
  fprintf(_file, "%d %d %lf %lf\n", nx, ny, dx, dy);
  fprintf(_file, "%lf %lf %lf %lf\n", xmin, xmax, ymin, ymax);

  int necr = 0;
  for (int iy = 0; iy < ny; iy++)
  {
    if (iy % _nsampley != 0) continue;
    for (int ix = 0; ix < nx; ix++)
    {
      if (ix % _nsamplex != 0) continue;
      indg[0] = ix;
      indg[1] = iy;
      int iech = _dbgrid->indiceToRank(indg);
      double value = _dbgrid->getArray(iech, _cols[0]);
      if (FFFF(value)) value = 9999990.;
      fprintf(_file, "%10.3lf ", value);
      necr++;
      if (necr == 6)
      {
        fprintf(_file, "\n");
        necr = 0;
      }
    }
  }
  if (necr > 0) fprintf(_file, "\n");

  _fileClose();
  return 0;
}
