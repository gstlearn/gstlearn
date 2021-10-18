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
#include "Basic/Table.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "Basic/AException.hpp"

Table::Table(int nrows, int ncols)
  : ASerializable(),
    _stats()
{
  init (nrows, ncols);
}

Table::Table(const VectorVectorDouble& table)
  : ASerializable(),
    _stats(table)
{
}


Table::Table(const String& neutralFileName, bool verbose)
    : ASerializable(),
      _stats()
{
  if (deSerialize(neutralFileName, verbose))
    my_throw("Problem reading the Neutral File");
}


Table::Table(const Table &m)
    : ASerializable(m),
      _stats(m._stats)
{

}

Table& Table::operator=(const Table &m)
{
  if (this != &m)
  {
    _stats = m._stats;
  }
  return *this;
}

Table::~Table()
{

}

int Table::getRowNumber() const
{
  if (isEmpty()) return 0;
  return static_cast<int>(_stats[0].size());
}

int Table::getColNumber() const
{
  if (isEmpty()) return 0;
  return static_cast<int>(_stats.size());
}

VectorDouble Table::getCol(int icol) const
{
  if (! _isColValid(icol)) return VectorDouble();
  return _stats[icol];
}

VectorDouble Table::getRow(int irow) const
{
  if (! _isRowValid(irow)) return VectorDouble();
  int ncols = getColNumber();
  VectorDouble vec(ncols);
  for (int icol = 0; icol < ncols; icol++)
    vec[icol] = _stats[icol][irow];
  return vec;
}

void Table::resize(int nrows, int ncols)
{
  int nrows_cur = 0;
  if (isEmpty())
    _stats.resize(ncols);
  else
    nrows_cur = getRowNumber();

  if (nrows >= nrows_cur)
  {
    int ncols = getColNumber();
    for (int icol = 0; icol < ncols; icol++)
      _stats[icol].resize(nrows);
  }
}

double Table::getValue(int irow, int icol) const
{
  if (icol < 0 || icol >= getColNumber()) return 0.;
  if (irow < 0 || irow >= getRowNumber()) return 0.;
  return _stats[icol][irow];
}

VectorDouble Table::getRange(int icol) const
{
  VectorDouble vec = getCol(icol);
  if (vec.empty()) return VectorDouble();
  VectorDouble limits(2);
  limits[0] = ut_vector_min(vec);
  limits[1] = ut_vector_max(vec);
  return limits;
}

VectorDouble Table::getRange() const
{
  int ncols = getColNumber();
  VectorDouble limits(2);
  limits[0] =  1.e30;
  limits[1] = -1.e30;
  for (int icol = 0; icol < ncols; icol++)
  {
    VectorDouble local = getRange(icol);
    if (local[0] < limits[0]) limits[0] = local[0];
    if (local[1] > limits[1]) limits[1] = local[1];
  }
  return limits;
}

int Table::serialize(const String& filename, bool verbose) const
{

  /* Opening the Data file */

  if (_fileOpen(filename, "Table", "w", verbose)) return 1;

  /* Writing the header */

  _recordWrite("%d", getColNumber());
  _recordWrite("#", "Number of Columns");
  _recordWrite("%d", getRowNumber());
  _recordWrite("#", "Number of Rows");

  /* Writing the tail of the file */

  for (int irow = 0; irow < getRowNumber(); irow++)
  {
    for (int icol = 0; icol < getColNumber(); icol++)
    {
      _recordWrite("%lf", _stats[icol][irow]);
    }
   _recordWrite("\n");
  }

  // Close the Neutral file
  _fileClose(verbose);

  return 0;
}

int Table::deSerialize(const String& filename, bool verbose)
{
  int ncols, nrows;
  double value;

  /* Opening the Data file */

  if (_fileOpen(filename, "Table", "r", verbose)) return 1;
  int error = 1;

  /* Decoding the header */

  if (_recordRead("Number of Columns", "%d", &ncols)) goto label_end;
  if (_recordRead("Number of Rows", "%d", &nrows)) goto label_end;

  _stats.clear();
  _stats.resize(ncols);

  /* Loop on the lines */

  for (int irow = 0; irow < nrows; irow++)
  {
    for (int icol = 0; icol < ncols; icol++)
    {
      if (_recordRead("Numerical value", "%lf", &value)) goto label_end;
      _stats[icol].push_back(value);
    }
  }

  error = 0;

  label_end:
  _fileClose(verbose);
  if (error) _stats.clear();
  return error;
}

/**
 * Print the contents of the statistics
 */
void Table::display(int isimu) const
{
  if (_stats.empty()) return;

  mestitle(1,"Statistics printout for Simulation %d",isimu);

  int ncols = getColNumber();
  int nrows = getRowNumber();

  for (int irow = 0; irow < nrows; irow++)
  {
    for (int icol = 0; icol < ncols; icol++)
      message(" %10.3lf",_stats[icol][irow]);
    message("\n");
  }
}

/**
 * Plot the contents of the statistics
 */
void Table::plot(int isimu) const
{
  if (_stats.empty()) return;
  String filename = incrementStringVersion("GibbsStats",isimu+1);
  serialize(filename,false);
}

bool Table::_isColValid(int icol) const
{
  int ncols = getColNumber();
  if (icol < 0 || icol >= ncols)
  {
    mesArg("Table Column", icol, ncols);
    return false;
  }
  return true;
}

bool Table::_isRowValid(int irow) const
{
  int nrows = getRowNumber();
  if (irow < 0 || irow >= nrows)
  {
    mesArg("Table Row", irow, nrows);
    return false;
  }
  return true;
}
