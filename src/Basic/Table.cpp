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
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"

Table::Table(int nrows, int ncols)
  : ASerializable(),
    AStringable(),
    _tab(),
    _rowNames(),
    _colNames()
{
  init(nrows, ncols);
}

Table::Table(const Table &m)
    : ASerializable(m),
      AStringable(m),
      _tab(m._tab),
      _rowNames(m._rowNames),
      _colNames(m._colNames)
{

}

Table& Table::operator=(const Table &m)
{
  if (this != &m)
  {
    ASerializable::operator=(m);
    AStringable::operator=(m);
    _tab = m._tab;
    _rowNames = m._rowNames;
    _colNames = m._colNames;
  }
  return *this;
}

Table::~Table()
{
}

int Table::resetFromArray(const VectorVectorDouble& table, bool flagByRow)
{
  _tab.reset(table, flagByRow);
  _rowNames.clear();
  _colNames.clear();
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

Table* Table::createFromArray(const VectorVectorDouble& tabin, bool flagByRow)
{
  Table* table = new Table();
  if (table->resetFromArray(tabin, flagByRow))
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
  return _tab.getNRows();
}

int Table::getColumnNumber() const
{
  if (isEmpty()) return 0;
  return _tab.getNCols();
}

/**
 * Extract a Column from a Table, stripping the TEST values
 * @param icol Rank of the target column
 * @return
 */
VectorDouble Table::getColumn(int icol) const
{
  if (! _isColValid(icol)) return VectorDouble();
  VectorDouble vec = _tab.getColumn(icol);
  return vec;
}

/**
 * Extract a Row, stripping the TEST values
 * @param irow Rank of the target row
 * @return
 */
VectorDouble Table::getRow(int irow) const
{
  if (! _isRowValid(irow)) return VectorDouble();
  return _tab.getRow(irow);
}

void Table::update(int irow, int icol, double value)
{
  _tab.setValue(irow, icol, value);
}

void Table::increment(int irow, int icol, double value)
{
  _tab.setValue(irow, icol, _tab.getValue(irow, icol) + value);
}

void Table::addRow()
{
  MatrixRectangular statsSave(_tab);
  int nrows = _tab.getNRows();
  int ncols = _tab.getNCols();

  _tab.reset(nrows+1, ncols);
  for (int irow=0; irow< nrows; irow++)
    for (int icol=0; icol<ncols; icol++)
    {
      _tab.setValue(irow, icol, statsSave.getValue(irow, icol));
    }
}

double Table::getValue(int irow, int icol) const
{
  if (icol < 0 || icol >= getColumnNumber()) return 0.;
  if (irow < 0 || irow >= getRowNumber()) return 0.;
  return _tab.getValue(irow, icol);
}

void Table::setValue(int irow, int icol, double value)
{
  if (icol < 0 || icol >= getColumnNumber()) return;
  if (irow < 0 || irow >= getRowNumber()) return;
  _tab.setValue(irow, icol, value);
}

VectorDouble Table::getRange(int icol) const
{
  VectorDouble vec = getColumn(icol);
  if (vec.empty()) return VectorDouble();
  VectorDouble limits(2);
  limits[0] = VH::minimum(vec);
  limits[1] = VH::maximum(vec);
  return limits;
}

VectorDouble Table::getAllRange() const
{
  int ncols = getColumnNumber();
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
  ret = ret && _recordWrite<int>(os, "Number of Columns", getColumnNumber());
  ret = ret && _recordWrite<int>(os, "Number of Rows", getRowNumber());

  /* Writing the tail of the file */

  for (int irow = 0; ret && irow < getRowNumber(); irow++)
  {
    for (int icol = 0; ret && icol < getColumnNumber(); icol++)
    {
      ret = ret && _recordWrite<double>(os, "", _tab.getValue(irow, icol));
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

  _tab.reset(nrows, ncols);

  /* Loop on the lines */

  for (int irow = 0; ret && irow < nrows; irow++)
  {
    for (int icol = 0; ret && icol < ncols; icol++)
    {
      ret = ret && _recordRead<double>(is, "Numerical value", value);
      if (ret) _tab.setValue(irow, icol, value);
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
  if (_tab.isEmpty()) return sstr.str();

  sstr << toTitle(1, "Table contents");
  int ncols = getColumnNumber();
  int nrows = getRowNumber();
  sstr << "- Number of Rows    = " << nrows << std::endl;
  sstr << "- Number of Columns = " << ncols << std::endl;
  sstr << std::endl;

  // Print optional header (using Column names if defined)

  if (! _colNames.empty())
  {
    if (! _rowNames.empty()) sstr << toStr(" ");
    for (int icol = 0; icol < ncols; icol++)
      sstr << " " << toStr(_colNames[icol]);
    sstr << std::endl;
  }

  // Print the contents of the table
  for (int irow = 0; irow < nrows; irow++)
  {
    if (! _rowNames.empty()) sstr << toStr(_rowNames[irow]);
    for (int icol = 0; icol < ncols; icol++)
    {
      sstr << " " << toDouble(_tab.getValue(irow, icol));
    }
    sstr << std::endl;
  }
  return sstr.str();
}

/**
 * Plot the contents of the statistics
 */
void Table::plot(int isimu) const
{
  if (_tab.isEmpty()) return;
  String filename = incrementStringVersion("TableStats",isimu+1);
  (void) dumpToNF(filename);
}

bool Table::_isColValid(int icol) const
{
  int ncols = getColumnNumber();
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

void Table::setColumnName(int icol, const String& name)
{
  if (! _isColValid(icol)) return;
  int ncols = getColumnNumber();
  if (_colNames.empty())
    _colNames.resize(ncols, "  ");
  _colNames[icol] = name;
}

void Table::setRowName(int irow, const String& name)
{
  if (! _isRowValid(irow)) return;
  int nrows = getRowNumber();
  if (_rowNames.empty())
    _rowNames.resize(nrows, "  ");
  _rowNames[irow] = name;
}

void Table::fill(double valinit)
{
  for (int irow = 0; irow < getRowNumber(); irow++)
    for (int icol = 0; icol < getColumnNumber(); icol++)
      _tab.setValue(irow, icol, valinit);
}
