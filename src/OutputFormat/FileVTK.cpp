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
#include "geoslib_f.h"
#include "geoslib_old_f.h"

#include "OutputFormat/FileVTK.hpp"
#include "OutputFormat/AOF.hpp"
#include "OutputFormat/vtk.h"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"

#include <string.h>

FileVTK::FileVTK(const char* filename, const Db* db)
  : AOF(filename, db)
  , _flagBinary(false)
  , _factx(1)
  , _facty(1)
  , _factz(1)
  , _factvar(1)
{
}

FileVTK::FileVTK(const FileVTK& r)
    : AOF(r),
      _flagBinary(r._flagBinary),
      _factx(r._factx),
      _facty(r._facty),
      _factz(r._factz),
      _factvar(r._factvar)
{
}

FileVTK& FileVTK::operator=(const FileVTK& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
    _flagBinary = r._flagBinary;
    _factx = r._factx;
    _facty = r._facty;
    _factz = r._factz;
    _factvar = r._factvar;
  }
  return *this;
}

FileVTK::~FileVTK()
{
}

int FileVTK::writeInFile()
{
  int dims[3];

  if (_fileWriteOpen()) return 1;

  /* Preliminary checks */

  int ndim = _db->getNDim();
  int ncol = (int) _cols.size();
  int nech = _db->getSampleNumber();
  int nactive = _db->getSampleNumber(true);
  bool flag_grid = _db->isGrid();

  /* Define the reading parameters */

  if (flag_grid)
    for (int idim = 0; idim < 3; idim++)
      dims[idim] = (idim < ndim) ? _dbgrid->getNX(idim) : 1;

  /* Core allocation */

  VectorInt vardim(ncol);
  VectorInt center(ncol);
  for (int icol = 0; icol < ncol; icol++)
  {
    vardim[icol] = 1;
    center[icol] = 1;
  }

  float** tab = (float**) mem_alloc(sizeof(float*) * ncol, 1);
  for (int icol = 0; icol < ncol; icol++)
  {
    if (flag_grid)
      tab[icol] = (float*) mem_alloc(sizeof(float) * nech, 1);
    else
      tab[icol] = (float*) mem_alloc(sizeof(float) * nactive, 1);
  }

  VectorFloat xcoor;
  VectorFloat ycoor;
  VectorFloat zcoor;
  VectorFloat points;
  if (flag_grid)
  {
    xcoor.resize(dims[0]);
    xcoor[0] = 0.;
    if (dims[0] > 1)
      for (int i = 0; i < dims[0]; i++)
        xcoor[i] = (float) (_factx * (_dbgrid->getX0(0) + i * _dbgrid->getDX(0)));
    ycoor.resize(dims[1]);
    ycoor[0] = 0.;
    if (dims[1] > 1)
      for (int i = 0; i < dims[1]; i++)
        ycoor[i] = (float) (_facty * (_dbgrid->getX0(1) + i * _dbgrid->getDX(1)));
    zcoor.resize(dims[2]);
    zcoor[0] = 0.;
    if (dims[2] > 1)
      for (int i = 0; i < dims[2]; i++)
        zcoor[i] = (float) (_factz * (_dbgrid->getX0(2) + i * _dbgrid->getDX(2)));
  }
  else
  {
    points.resize(3 * nactive);
  }

  /* Read the coordinates (for points only) */

  if (! flag_grid)
  {
    int ecr = 0;
    for (int iech = 0; iech < nech; iech++)
    {
      if (! _db->isActive(iech)) continue;
      for (int idim = 0; idim < 3; idim++)
      {
        int fact = 1;
        if (idim == 0) fact = _factx;
        if (idim == 1) fact = _facty;
        if (idim == 2) fact = _factz;
        points[ecr++] = (idim < ndim) ? (float) (fact * _db->getCoordinate(iech, idim)) : 0.;
      }
    }
  }

  /* Load the array */

  for (int icol = 0; icol < ncol; icol++)
  {
    if (!flag_grid)
    {
      int ecr = 0;
      for (int iech = 0; iech < nech; iech++)
        if (_db->isActive(iech))
        {
          double value = (float) (_db->getArray(iech, _cols[icol]));
          if (FFFF(value))
            tab[icol][ecr] = (float) (TEST);
          else
            tab[icol][ecr] = (float) (_factvar * value);
          ecr++;
        }
    }
    else
    {
      int ecr = 0;
      for (int iz = 0; iz < dims[2]; iz++)
        for (int iy = 0; iy < dims[1]; iy++)
          for (int ix = 0; ix < dims[0]; ix++)
          {
            int iad = ix + dims[0] * (iy + dims[1] * iz);
            if (_dbgrid->isActive(iad))
            {
              double value = (float) (_dbgrid->getValueByColIdx(iad, _cols[icol]));
              if (FFFF(value))
                tab[icol][ecr] = (float) (TEST);
              else
                tab[icol][ecr] = (float) (_factvar * value);
            }
            else
              tab[icol][ecr] = (float) (TEST);
            ecr++;
          }
    }
  }

  VectorString names = _db->getNamesByColIdx(_cols);
  std::vector<char*> vc = util_vs_to_vs(names);

  /* Write the file */

  if (flag_grid)
    write_rectilinear_mesh(getFilename().c_str(), _flagBinary,
                           dims, xcoor.data(), ycoor.data(), zcoor.data(), ncol,
                           vardim.data(), center.data(), vc.data(), tab);
  else
    write_point_mesh(getFilename().c_str(), _flagBinary,
                     nactive, points.data(), ncol, vardim.data(), vc.data(), tab);

  for (int icol = 0; icol < ncol; icol++)
    tab[icol] = (float*) mem_free((char* ) tab[icol]);
  tab = (float**) mem_free((char* ) tab);

  _fileClose();
  return 0;
}
