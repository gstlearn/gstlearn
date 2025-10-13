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
#include "Matrix/Table.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/String.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
Table::Table(Id nrow, Id ncol, bool skip_title, bool skip_description)
  : MatrixDense(nrow, ncol)
  , ASerializable()
  , _title()
  , _rowNames()
  , _colNames()
  , _skipTitle(skip_title)
  , _skipDescription(skip_description)
{
  reset(nrow, ncol);
}

Table::Table(const Table& m)
  : MatrixDense(m)
  , ASerializable(m)
  , _title(m._title)
  , _rowNames(m._rowNames)
  , _colNames(m._colNames)
  , _skipTitle(m._skipTitle)
  , _skipDescription(m._skipDescription)
{
}

Table& Table::operator=(const Table& m)
{
  if (this != &m)
  {
    MatrixDense::operator=(m);
    ASerializable::operator=(m);
    _title           = m._title;
    _rowNames        = m._rowNames;
    _colNames        = m._colNames;
    _skipTitle       = m._skipTitle;
    _skipDescription = m._skipDescription;
  }
  return *this;
}

Table::~Table()
{
}

void Table::reset(Id nrows, Id ncols)
{
  MatrixDense::reset(nrows, ncols);
  _clearDecoration();
}

void Table::_clearDecoration()
{
  _rowNames.clear();
  _colNames.clear();
}

Table* Table::create(Id nrow, Id ncol)
{
  return new Table(nrow, ncol);
}

Table* Table::createFromNames(const VectorString& rownames,
                              const VectorString& colnames)
{
  Id nrow    = static_cast<Id>(rownames.size());
  Id ncol    = static_cast<Id>(colnames.size());
  auto* table = new Table(nrow, ncol);
  table->setRowNames(rownames);
  table->setColumnNames(colnames);
  return table;
}

Table* Table::createFromNF(const String& NFFilename, bool verbose)
{
  auto* table = new Table();
  if (table->_fileOpenAndDeserialize(NFFilename, verbose)) return table;
  delete table;
  return nullptr;
}

Table* Table::createFromTable(const Table& table)
{
  return new Table(table);
}

VectorDouble Table::getRange(Id icol) const
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
  auto ncols = getNCols();
  VectorDouble limits(2);
  limits[0] = MAXIMUM_BIG;
  limits[1] = MINIMUM_BIG;
  for (Id icol = 0; icol < ncols; icol++)
  {
    VectorDouble local = getRange(icol);
    if (local[0] < limits[0]) limits[0] = local[0];
    if (local[1] > limits[1]) limits[1] = local[1];
  }
  return limits;
}

bool Table::_serializeAscii(std::ostream& os, bool /*verbose*/) const
{
  bool ret = true;
  ret      = ret && _recordWrite<Id>(os, "Number of Columns", getNCols());
  ret      = ret && _recordWrite<Id>(os, "Number of Rows", getNRows());

  /* Writing the tail of the file */

  for (Id irow = 0; ret && irow < getNRows(); irow++)
  {
    for (Id icol = 0; ret && icol < getNCols(); icol++)
    {
      ret = ret && _recordWrite<double>(os, "", getValue(irow, icol));
    }
    ret = ret && _commentWrite(os, "");
  }
  return ret;
}

bool Table::_deserializeAscii(std::istream& is, bool /*verbose*/)
{
  Id nrows    = 0;
  Id ncols    = 0;
  double value = 0.;

  bool ret = true;
  ret      = ret && _recordRead<Id>(is, "Number of Columns", ncols);
  ret      = ret && _recordRead<Id>(is, "Number of Rows", nrows);
  if (!ret) return false;

  reset(nrows, ncols);

  /* Loop on the lines */

  for (Id irow = 0; ret && irow < nrows; irow++)
  {
    for (Id icol = 0; ret && icol < ncols; icol++)
    {
      ret = ret && _recordRead<double>(is, "Numerical value", value);
      if (ret) setValue(irow, icol, value);
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
  if (empty()) return sstr.str();
  auto ncols = getNCols();
  auto nrows = getNRows();

  // Title
  if (!_skipTitle)
  {
    if (!_title.empty())
      sstr << toTitle(1, _title.c_str());
  }

  // Description
  if (!_skipDescription)
  {
    sstr << "- Number of Rows    = " << nrows << std::endl;
    sstr << "- Number of Columns = " << ncols << std::endl;
    sstr << std::endl;
  }

  // For displaying the Row names, find the optimal dimension
  Id rowLengthMax = 1;
  if (!_rowNames.empty())
  {
    for (Id irow = 0; irow < nrows; irow++)
    {
      Id rowLength = static_cast<Id>(_rowNames[irow].size());
      if (rowLength > rowLengthMax) rowLengthMax = rowLength;
    }
  }

  // Print optional header (using Column names if defined)
  if (!_colNames.empty())
  {
    if (!_rowNames.empty()) sstr << toStr(" ", EJustify::fromKey("RIGHT"), rowLengthMax);
    for (Id icol = 0; icol < ncols; icol++)
      sstr << " " << toStr(_colNames[icol]);
    sstr << std::endl;
  }

  // Print the contents of the table
  for (Id irow = 0; irow < nrows; irow++)
  {
    if (!_rowNames.empty()) sstr << toStr(_rowNames[irow], EJustify::fromKey("RIGHT"), rowLengthMax);
    for (Id icol = 0; icol < ncols; icol++)
    {
      sstr << " " << toDouble(getValue(irow, icol));
    }
    sstr << std::endl;
  }
  return sstr.str();
}

/**
 * Plot the contents of the statistics
 */
void Table::plot(Id isimu) const
{
  if (empty()) return;
  String filename = incrementStringVersion("TableStats", isimu + 1);
  (void)dumpToNF(filename);
}

void Table::setColumnNames(const VectorString& colNames)
{
  if (getNCols() != static_cast<Id>(colNames.size()))
  {
    messerr("The size of 'colNames' (%d) does not match the number of columns (%d)",
            static_cast<Id>(colNames.size()), getNCols());
    return;
  }
  _colNames = colNames;
}

void Table::setColumnName(Id icol, const String& name)
{
  if (!_isColumnValid(icol)) return;
  auto ncols = getNCols();
  if (_colNames.empty())
    _colNames.resize(ncols, "  ");
  _colNames[icol] = name;
}

void Table::setRowName(Id irow, const String& name)
{
  if (!_isRowValid(irow)) return;
  auto nrows = getNRows();
  if (_rowNames.empty())
    _rowNames.resize(nrows, "  ");
  _rowNames[irow] = name;
}

void Table::setRowNames(const VectorString& rowNames)
{
  if (getNRows() != static_cast<Id>(rowNames.size()))
  {
    messerr("The size of 'rowNames' (%d) does not match the number of rows (%d)",
            static_cast<Id>(rowNames.size()), getNRows());
    return;
  }
  _rowNames = rowNames;
}

String Table::getColumnName(Id icol) const
{
  if (!_isColumnValid(icol)) return String();
  return _colNames[icol];
}

String Table::getRowName(Id irow) const
{
  if (!_isRowValid(irow)) return String();
  return _rowNames[irow];
}
#ifdef HDF5
bool Table::_deserializeH5(H5::Group& grp, [[maybe_unused]] bool verbose)
{
  auto tableG = SerializeHDF5::getGroup(grp, "Table");
  if (!tableG)
  {
    return false;
  }

  /* Read the grid characteristics */
  bool ret = true;

  Id ncols = 0;
  Id nrows = 0;
  VectorDouble values;

  ret = ret && SerializeHDF5::readValue(*tableG, "NCols", ncols);
  ret = ret && SerializeHDF5::readValue(*tableG, "NRows", nrows);
  ret = ret && SerializeHDF5::readVec(*tableG, "Values", values);
  if (ret)
  {
    reset(nrows, ncols);
    setValues(values);
  }

  return ret;
}

bool Table::_serializeH5(H5::Group& grp, [[maybe_unused]] bool verbose) const
{
  auto tableG = grp.createGroup("Table");

  bool ret = true;

  ret = ret && SerializeHDF5::writeValue(tableG, "NCols", getNCols());
  ret = ret && SerializeHDF5::writeValue(tableG, "NRows", getNRows());
  ret = ret && SerializeHDF5::writeVec(tableG, "Values", getValues());

  return ret;
}
#endif
} // namespace gstlrn
