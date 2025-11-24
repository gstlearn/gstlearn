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
#include "Basic/StringHelper.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/String.hpp"
#include "Enum/ECst.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <locale>
#include <regex>
#include <sstream>

#define CASE_DOUBLE 0
#define CASE_REAL   1
#define CASE_INT    2
#define CASE_COL    3
#define CASE_ROW    4

namespace gstlrn
{
Id _getColumnRank()
{
  return static_cast<Id>(OptCst::query(ECst::NTRANK));
}
Id _getColumnName()
{
  return static_cast<Id>(OptCst::query(ECst::NTNAME));
}
Id _getColumnSize(Id localSize, Id nColumns)
{
  if (localSize > 0) return localSize;
  Id size = 1 + static_cast<Id>(OptCst::query(ECst::NTCAR));
  if (nColumns == 1) return size;
  return size * nColumns; // account for spaces between columns
}
Id _getMaxNCols()
{
  return static_cast<Id>(OptCst::query(ECst::NTCOL));
}
Id _getMaxNRows()
{
  return static_cast<Id>(OptCst::query(ECst::NTROW));
}
Id _getNBatch()
{
  return static_cast<Id>(OptCst::query(ECst::NTBATCH));
}
Id _getDecimalNumber()
{
  return static_cast<Id>(OptCst::query(ECst::NTDEC));
}
double _getThresh()
{
  Id ndec       = static_cast<Id>(OptCst::query(ECst::NTDEC));
  double thresh = (0.5 * pow(10, -ndec));
  return thresh;
}

String _toStrRowColumn(Id icase, Id value, Id flagAdd)
{
  std::stringstream sstr;
  I32 rank  = static_cast<I32>(_getColumnRank());
  I32 width = static_cast<I32>(_getColumnSize() - _getColumnRank() - 1);
  sstr << std::setw(width) << std::right;
  if (icase == CASE_ROW)
  {
    if (!flagAdd)
      sstr << "[" << std::setw(rank) << value << ",]";
    else
      sstr << "[" << std::setw(rank) << value << "+]";
  }
  else
  {
    if (!flagAdd)
      sstr << "[," << std::setw(rank) << value << "]";
    else
      sstr << "[ " << std::setw(rank) << value << "]";
  }
  return sstr.str();
}

// Functions for beautifying a suite of values
String _toStrRowHeader(const VectorString& rownames, Id iy, Id rowSize)
{
  std::stringstream sstr;
  if (!rownames.empty())
    sstr << toStr(rownames[iy], -1, rowSize);
  else
    sstr << _toStrRowColumn(CASE_ROW, iy, false);
  return sstr.str();
}

String _toStrColumnHeaders(const VectorString& colnames, Id colfrom, Id colto, Id colSize)
{
  std::stringstream sstr;
  if (!colnames.empty())
  {
    // By Names
    sstr << toStr(" ", 1) << " ";
    for (Id ix = colfrom; ix < colto; ix++)
      sstr << toStr(colnames[ix], 1, colSize);
    sstr << std::endl;
  }
  else
  {
    // By Numbers
    sstr << toStr(" ", 1) << " ";
    for (Id ix = colfrom; ix < colto; ix++)
      sstr << _toStrRowColumn(CASE_COL, ix, false);
    sstr << std::endl;
  }
  return sstr.str();
}

String _toStrTrailer(Id ncols, Id nrows, Id ncols_util, Id nrows_util)
{
  std::stringstream sstr;

  bool one_used = (ncols != ncols_util || nrows != nrows_util);
  bool all_used = (ncols != ncols_util && nrows != nrows_util);

  if (one_used) sstr << "(";

  if (ncols != ncols_util)
  {
    if (ncols == ncols_util)
      sstr << "Ncols=" << ncols;
    else
      sstr << "Ncols=" << ncols_util << "[from " << ncols << "]";
  }

  if (all_used) sstr << ",";

  if (nrows != nrows_util)
  {
    if (nrows == nrows_util)
      sstr << "Nrows=" << nrows;
    else
      sstr << "Nrows=" << nrows_util << "[from " << nrows << "]";
  }

  if (one_used) sstr << ")" << std::endl;
  return sstr.str();
}
} // namespace gstlrn