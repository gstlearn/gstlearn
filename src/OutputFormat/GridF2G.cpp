/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"

#include "OutputFormat/AOF.hpp"
#include "OutputFormat/GridF2G.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>

#define F2G(ix,iy,iz,icol)  (tab[(ix) + nx[0] * ((iy) + nx[1] * ((iz) + nx[2] * (icol)))])

GridF2G::GridF2G(const char* filename, const Db* db)
  : AOF(filename, db)
{
}

GridF2G::GridF2G(const GridF2G& r)
    : AOF(r)
{
}

GridF2G& GridF2G::operator=(const GridF2G& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
  }
  return *this;
}

GridF2G::~GridF2G()
{
}

DbGrid* GridF2G::readGridFromFile()
{
  DbGrid* dbgrid = nullptr;
  char string[100], refchar[100], valtest[10], valread[10];
  VectorInt nx(3);
  VectorDouble x0(3);
  VectorDouble dx(3);
  VectorDouble angles(3);
  VectorString names;
  int ndim, version, ncol;
  double dum, value;

  /* Open the file */

  if (_fileReadOpen()) return dbgrid;

  /* Initializations */

  for (int idim = 0; idim < 3; idim++)
  {
    nx[idim] = 1;
    x0[idim] = 0.;
    dx[idim] = 1.;
    angles[idim] = 0.;
  }

  /* Read the header */

  if (_record_read(_file, "%s", string)) return dbgrid;
  (void) gslStrcpy(refchar, "F2G_DIM");
  if (strcmp(string, refchar)) return dbgrid;
  if (_record_read(_file, "%d", &ndim)) return dbgrid;

  if (_record_read(_file, "%s", string)) return dbgrid;
  (void) gslStrcpy(refchar, "F2G_VERSION");
  if (strcmp(string, refchar)) return dbgrid;
  if (_record_read(_file, "%d", &version)) return dbgrid;

  if (_record_read(_file, "%s", string)) return dbgrid;
  (void) gslStrcpy(refchar, "F2G_LOCATION");
  if (strcmp(string, refchar)) return dbgrid;
  for (int idim = 0; idim < 3; idim++) // Always three parameters
    if (_record_read(_file, "%lf", &x0[idim])) return dbgrid;

  if (_record_read(_file, "%s", string)) return dbgrid;
  (void) gslStrcpy(refchar, "F2G_ROTATION");
  if (strcmp(string, refchar)) return dbgrid;
  if (_record_read(_file, "%lf", angles[0])) return dbgrid;

  if (_record_read(_file, "%s", string)) return dbgrid;
  (void) gslStrcpy(refchar, "F2G_ORIGIN");
  if (strcmp(string, refchar)) return dbgrid;
  for (int idim = 0; idim < ndim; idim++)
    if (_record_read(_file, "%lf", &dum)) return dbgrid;

  if (_record_read(_file, "%s", string)) return dbgrid;
  (void) gslStrcpy(refchar, "F2G_NB_NODES");
  if (strcmp(string, refchar)) return dbgrid;
  for (int idim = 0; idim < ndim; idim++)
    if (_record_read(_file, "%d", &nx[idim])) return dbgrid;

  if (_record_read(_file, "%s", string)) return dbgrid;
  (void) gslStrcpy(refchar, "F2G_LAGS");
  if (strcmp(string, refchar)) return dbgrid;
  for (int idim = 0; idim < ndim; idim++)
    if (_record_read(_file, "%lf", &dx[idim])) return dbgrid;

  if (_record_read(_file, "%s", string)) return dbgrid;
  (void) gslStrcpy(refchar, "F2G_ORDER");
  if (strcmp(string, refchar)) return dbgrid;

  // We need to read the three next strings (orders)
  // Only the order +Y +X +Z is interfaced
  if (_record_read(_file, "%s", string)) return dbgrid;
  if (strcmp(string, "+Y")) return dbgrid;
  if (_record_read(_file, "%s", string)) return dbgrid;
  if (strcmp(string, "+X")) return dbgrid;
  if (_record_read(_file, "%s", string)) return dbgrid;
  if (strcmp(string, "+Z")) return dbgrid;

  if (_record_read(_file, "%s", string)) return dbgrid;
  (void) gslStrcpy(refchar, "F2G_NB_VARIABLES");
  if (strcmp(string, refchar)) return dbgrid;
  if (_record_read(_file, "%d", &ncol)) return dbgrid;

  for (int i = 0; i < ncol; i++)
  {
    if (_record_read(_file, "%s", string)) return dbgrid;
    // Variable Name
    (void) gslSPrintf(refchar, "F2G_VARIABLE_%d", i + 1);
    if (strcmp(string, refchar)) return dbgrid;
    // We need to read the Name even if ignored
    if (_record_read(_file, "%s", string)) return dbgrid;
    names.push_back(string);
    if (_record_read(_file, "%s", string)) return dbgrid;
    // NA value
    (void) gslSPrintf(refchar, "F2G_UNDEFINED_%d", i + 1);
    if (strcmp(string, refchar)) return dbgrid;
    if (_record_read(_file, "%s", valtest)) return dbgrid;
  }

  if (_record_read(_file, "%s", string)) return dbgrid;
  (void) gslStrcpy(refchar, "F2G_VALUES");
  if (strcmp(string, refchar)) return dbgrid;

  int size = nx[0] * nx[1] * nx[2];
  VectorDouble tab(size * ncol);
  for (int i = 0; i < size * ncol; i++) tab[i] = TEST;

  int nused = 0;
  for (int iz = 0; iz < nx[2]; iz++)
    for (int ix = 0; ix < nx[0]; ix++)
      for (int iy = 0; iy < nx[1]; iy++)
        for (int icol = 0; icol < ncol; icol++)
        {
          if (_record_read(_file, "%s", valread)) return dbgrid;
          if (!strcmp(valread, valtest))
            value = TEST;
          else
          {
            value = atof(valread);
            nused++;
          }
          F2G(ix,iy,iz,icol) = value;
        }

  /* Patch the origin of the grid along vertical (to match RGeostats) */

  x0[2] = x0[2] - dx[2] / 2.;

  dbgrid = new DbGrid();
  dbgrid->reset(nx,dx,x0,angles,ELoadBy::SAMPLE,tab,names);

  // Close the file

  _fileClose();

  return dbgrid;
}
