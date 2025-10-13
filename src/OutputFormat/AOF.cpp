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
#include "OutputFormat/AOF.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Db/DbGrid.hpp"

#include <cstdio>

namespace gstlrn
{
AOF::AOF(const String& filename, const Db* db)
  : _filename(filename)
  , _db(db)
  , _dbgrid(nullptr)
  , _cols()
  , _file(nullptr)
{
  if (db != nullptr) _dbgrid = dynamic_cast<const DbGrid*>(db);
}

AOF::AOF(const AOF& r)
  : _filename(r._filename)
  , _db(r._db)
  , _dbgrid(r._dbgrid)
  , _cols(r._cols)
  , _file(r._file)
{
}

AOF& AOF::operator=(const AOF& r)
{
  if (this != &r)
  {
    _filename = r._filename;
    _db       = r._db;
    _dbgrid   = r._dbgrid;
    _cols     = r._cols;
    _file     = r._file;
  }
  return *this;
}

AOF::~AOF()
{
}

bool AOF::isValidForGrid() const
{
  if (!mustBeGrid()) return true;
  if (_dbgrid == nullptr)
  {
    messerr("This function requires a Grid organization");
    return false;
  }
  return true;
}

bool AOF::isValidForVariable() const
{
  Id ncol = static_cast<Id>(_cols.size());
  if (mustBeOneVariable() && ncol > 1)
  {
    messerr("This function requires a single Variable but ncol = %d", ncol);
    return false;
  }
  return true;
}

bool AOF::isValidForNDim() const
{
  Id ndim = _dbgrid->getNDim();
  if (!mustBeForNDim(ndim))
  {
    messerr("This function is not valid for the Space Dimension (%d)", ndim);
    return false;
  }
  return true;
}

bool AOF::isValidForRotation() const
{
  Id ndim = _dbgrid->getNDim();

  Id mode = 0;
  if (_dbgrid->isGridRotated())
  {
    mode                = 1;
    VectorDouble angles = _dbgrid->getAngles();
    for (Id idim = 1; idim < ndim; idim++)
      if (ABS(angles[idim]) > 1.e-6) mode = 2;
  }
  if (!mustBeForRotation(mode))
  {
    messerr("This function is not compatible with Grid Rotation (mode=%d)", mode);
    return false;
  }
  return true;
}

Id AOF::_fileWriteOpen()
{
  _file = gslFopen(_filename.c_str(), "w");
  if (_file == nullptr)
  {
    messerr("Error when opening the file %s for writing", _filename.c_str());
    return (1);
  }
  return 0;
}

Id AOF::_fileReadOpen()
{
  _file = gslFopen(_filename.c_str(), "r");
  if (_file == nullptr)
  {
    messerr("Error when opening the file %s for reading", _filename.c_str());
    return (1);
  }
  return 0;
}

void AOF::_fileClose()
{
  if (_file != nullptr)
    fclose(_file);
  _file = nullptr;
}

void AOF::setCols(Id ncol, const Id* icols)
{
  _cols = VectorInt(ncol);
  for (Id icol = 0; icol < ncol; icol++)
    _cols[icol] = icols[icol];
}

void AOF::setCol(Id icol)
{
  _cols    = VectorInt(1);
  _cols[0] = icol;
}

bool AOF::isAuthorized() const
{
  if (_db == nullptr)
  {
    messerr("The argument 'db' must be provided");
    return false;
  }
  if (!isValidForGrid()) return false;
  if (!isValidForVariable()) return false;
  if (!isValidForNDim()) return false;
  if (!isValidForRotation()) return false;
  return true;
}
} // namespace gstlrn