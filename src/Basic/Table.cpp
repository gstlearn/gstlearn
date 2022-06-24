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
#include "Basic/Vector.hpp"
#include "Basic/AException.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"

Table::Table(int nrows, int ncols)
  : ASerializable(),
    AStringable(),
    _stats()
{
  init(nrows, ncols,true);
}

Table::Table(const Table &m)
    : ASerializable(m),
      AStringable(m),
      _stats(m._stats)
{

}

Table& Table::operator=(const Table &m)
{
  if (this != &m)
  {
    ASerializable::operator=(m);
    AStringable::operator=(m);
    _stats = m._stats;
  }
  return *this;
}

Table::~Table()
{
}

int Table::resetFromArray(const VectorVectorDouble& table)
{
  _stats = table;
  return 0;
}

Table* Table::create(int nrows, int ncols)
{
  return new Table(nrows, ncols);
}

Table* Table::createFromNF(const String& neutralFilename, bool verbose)
{
  Table* table = nullptr;
  std::ifstream is;
  table = new Table();
  bool success = false;
  if (table->_fileOpenRead(neutralFilename, is, verbose))
  {
    success =  table->deserialize(is, verbose);
  }
  if (! success)
  {
    delete table;
    table = nullptr;
  }
  return table;
}

Table* Table::createFromArray(const VectorVectorDouble& tabin)
{
  Table* table = new Table();
  if (table->resetFromArray(tabin))
  {
    messerr("Problem when loading a Table from Array");
    delete table;
    table = nullptr;
  }
  return table;
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

void Table::resize(int nrows, int ncols, bool zero)
{
  int nrows_cur = 0;
  if (isEmpty())
    _stats.resize(ncols);
  else
    nrows_cur = getRowNumber();

  if (nrows >= nrows_cur)
  {
    for (int icol = 0; icol < ncols; icol++)
    {
      if (zero)
        _stats[icol].resize(nrows,0.);
      else
        _stats[icol].resize(nrows);
    }
  }
}

double Table::getValue(int irow, int icol) const
{
  if (icol < 0 || icol >= getColNumber()) return 0.;
  if (irow < 0 || irow >= getRowNumber()) return 0.;
  return _stats[icol][irow];
}

void Table::setValue(int irow, int icol, double value)
{
  if (icol < 0 || icol >= getColNumber()) return;
  if (irow < 0 || irow >= getRowNumber()) return;
  _stats[icol][irow] = value;
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

VectorDouble Table::getAllRange() const
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

bool Table::_serialize(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret = ret && _recordWrite<int>(os, "Number of Columns", getColNumber());
  ret = ret && _recordWrite<int>(os, "Number of Rows", getRowNumber());

  /* Writing the tail of the file */

  for (int irow = 0; ret && irow < getRowNumber(); irow++)
  {
    for (int icol = 0; ret && icol < getColNumber(); icol++)
    {
      ret = ret && _recordWrite<double>(os, "", _stats[icol][irow]);
    }
    ret = ret && _commentWrite(os, "");
  }
  return ret;
}

bool Table::_deserialize(std::istream& is, bool /*verbose*/)
{
  int nrows = 0;
  int ncols = 0;
  double value = 0.;

  bool ret = true;
  ret = ret && _recordRead<int>(is, "Number of Columns", ncols);
  ret = ret && _recordRead<int>(is, "Number of Rows", nrows);
  if (! ret) return false;

  _stats.clear();
  _stats.resize(ncols);

  /* Loop on the lines */

  for (int irow = 0; ret && irow < nrows; irow++)
  {
    for (int icol = 0; ret && icol < ncols; icol++)
    {
      ret = ret && _recordRead<double>(is, "Numerical value", value);
      if (ret) _stats[icol].push_back(value);
    }
  }
  return ret;
}


/**
 * Print the contents of the statistics
 */
String Table::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (_stats.empty()) return sstr.str();

  sstr << toTitle(1, "Table contents");
  int ncols = getColNumber();
  int nrows = getRowNumber();
  sstr << "- Number of Rows    = " << nrows << std::endl;
  sstr << "- Number of Columns = " << ncols << std::endl;
  sstr << std::endl;

  for (int irow = 0; irow < nrows; irow++)
  {
    for (int icol = 0; icol < ncols; icol++)
      sstr << " " << toDouble(_stats[icol][irow]);
    sstr << std::endl;
  }
  return sstr.str();
}

/**
 * Plot the contents of the statistics
 */
void Table::plot(int isimu) const
{
  if (_stats.empty()) return;
  String filename = incrementStringVersion("TableStats",isimu+1);
  (void) dumpToNF(filename,false);
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
