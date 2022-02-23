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
#include "OutputFormat/AOF.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Basic/Vector.hpp"
#include "Db/DbGrid.hpp"
#include <stdio.h>

AOF::AOF(const char* filename, const Db* db)
  : _filename(filename)
  , _db(db)
  , _dbgrid(nullptr)
  , _cols()
  , _file(nullptr)
{
  _dbgrid = dynamic_cast<const DbGrid*>(db);
}

AOF::AOF(const AOF& r)
    : _filename(r._filename),
      _db(r._db),
      _dbgrid(r._dbgrid),
      _cols(r._cols),
      _file(r._file)
{
}

AOF& AOF::operator=(const AOF& r)
{
  if (this != &r)
  {
    _filename = r._filename;
    _db = r._db;
    _dbgrid = r._dbgrid;
    _cols = r._cols;
    _file = r._file;
  }
  return *this;
}

AOF::~AOF()
{
}

bool AOF::isValidForGrid() const
{
  if (! mustBeGrid()) return true;
  if (_dbgrid == nullptr)
  {
    messerr("This function requires a Grid organization");
    return false;
  }
  return true;
}

bool AOF::isValidForVariable() const
{
  int ncol = _cols.size();
  if (mustBeOneVariable() && ncol > 1)
  {
    messerr("This function requires a single Variable but ncol = %d",ncol);
    return false;
  }
  return true;
}

int AOF::_fileOpen()
{
  _file = gslFopen(_filename, "w");
  if (_file == nullptr)
  {
    messerr("Error when opening the file %s for writing", _filename);
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

void AOF::setCols(int ncol, int* icols)
{
  _cols = VectorInt(ncol);
  for (int icol = 0; icol < ncol; icol++)
    _cols[icol] = icols[icol];
}

void AOF::setCol(int icol)
{
  _cols = VectorInt(1);
  _cols[0] = icol;
}
