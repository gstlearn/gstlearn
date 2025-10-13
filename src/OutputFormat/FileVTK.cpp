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
#include "OutputFormat/FileVTK.hpp"
#include "OutputFormat/AOF.hpp"
#include "OutputFormat/vtk.h"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"

namespace gstlrn
{
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
  : AOF(r)
  , _flagBinary(r._flagBinary)
  , _factx(r._factx)
  , _facty(r._facty)
  , _factz(r._factz)
  , _factvar(r._factvar)
{
}

FileVTK& FileVTK::operator=(const FileVTK& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
    _flagBinary = r._flagBinary;
    _factx      = r._factx;
    _facty      = r._facty;
    _factz      = r._factz;
    _factvar    = r._factvar;
  }
  return *this;
}

FileVTK::~FileVTK()
{
}

Id FileVTK::writeInFile()
{
  Id dims[3];

  if (_fileWriteOpen()) return 1;

  /* Preliminary checks */

  Id ndim       = _db->getNDim();
  Id ncol       = static_cast<Id>(_cols.size());
  Id nech       = _db->getNSample();
  Id nactive    = _db->getNSample(true);
  bool flag_grid = _db->isGrid();

  /* Define the reading parameters */

  if (flag_grid)
    for (Id idim = 0; idim < 3; idim++)
      dims[idim] = (idim < ndim) ? _dbgrid->getNX(idim) : 1;

  /* Core allocation */

  VectorInt vardim(ncol);
  VectorInt center(ncol);
  for (Id icol = 0; icol < ncol; icol++)
  {
    vardim[icol] = 1;
    center[icol] = 1;
  }

  VectorVectorFloat tab(ncol);
  for (Id icol = 0; icol < ncol; icol++)
  {
    if (flag_grid)
      tab[icol].resize(nech);
    else
      tab[icol].resize(nactive);
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
      for (Id i = 0; i < dims[0]; i++)
        xcoor[i] = static_cast<float>(_factx * (_dbgrid->getX0(0) + i * _dbgrid->getDX(0)));
    ycoor.resize(dims[1]);
    ycoor[0] = 0.;
    if (dims[1] > 1)
      for (Id i = 0; i < dims[1]; i++)
        ycoor[i] = static_cast<float>(_facty * (_dbgrid->getX0(1) + i * _dbgrid->getDX(1)));
    zcoor.resize(dims[2]);
    zcoor[0] = 0.;
    if (dims[2] > 1)
      for (Id i = 0; i < dims[2]; i++)
        zcoor[i] = static_cast<float>(_factz * (_dbgrid->getX0(2) + i * _dbgrid->getDX(2)));
  }
  else
  {
    points.resize(3 * nactive);
  }

  /* Read the coordinates (for points only) */

  if (!flag_grid)
  {
    Id ecr = 0;
    for (Id iech = 0; iech < nech; iech++)
    {
      if (!_db->isActive(iech)) continue;
      for (Id idim = 0; idim < 3; idim++)
      {
        Id fact = 1;
        if (idim == 0) fact = _factx;
        if (idim == 1) fact = _facty;
        if (idim == 2) fact = _factz;
        points[ecr++] = (idim < ndim) ? static_cast<float>(fact * _db->getCoordinate(iech, idim)) : 0.;
      }
    }
  }

  /* Load the array */

  for (Id icol = 0; icol < ncol; icol++)
  {
    if (!flag_grid)
    {
      Id ecr = 0;
      for (Id iech = 0; iech < nech; iech++)
        if (_db->isActive(iech))
        {
          double value = static_cast<float>(_db->getArray(iech, _cols[icol]));
          if (FFFF(value))
            tab[icol][ecr] = static_cast<float>(TEST);
          else
            tab[icol][ecr] = static_cast<float>(_factvar * value);
          ecr++;
        }
    }
    else
    {
      Id ecr = 0;
      for (Id iz = 0; iz < dims[2]; iz++)
        for (Id iy = 0; iy < dims[1]; iy++)
          for (Id ix = 0; ix < dims[0]; ix++)
          {
            Id iad = ix + dims[0] * (iy + dims[1] * iz);
            if (_dbgrid->isActive(iad))
            {
              double value = static_cast<float>(_dbgrid->getValueByColIdx(iad, _cols[icol]));
              if (FFFF(value))
                tab[icol][ecr] = static_cast<float>(TEST);
              else
                tab[icol][ecr] = static_cast<float>(_factvar * value);
            }
            else
              tab[icol][ecr] = static_cast<float>(TEST);
            ecr++;
          }
    }
  }

  VectorString names = _db->getNamesByColIdx(_cols);
  std::vector<const char*> vc(names.size(), nullptr);
  for (Id i = 0; i < static_cast<Id>(names.size()); i++)
  {
    vc[i] = names[i].c_str();
  }

  /* Write the file */

  if (flag_grid)
    write_rectilinear_mesh(getFilename().c_str(), _flagBinary, dims,
                           xcoor.data(), ycoor.data(), zcoor.data(), ncol,
                           vardim.data(), center.data(), vc.data(), tab);
  else
    write_point_mesh(getFilename().c_str(), _flagBinary, nactive, points.data(),
                     ncol, vardim.data(), vc.data(), tab);

  _fileClose();
  return 0;
}
}
