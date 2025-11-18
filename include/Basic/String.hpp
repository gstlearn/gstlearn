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
#pragma once

#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

// TODO : add Namespace
#define SPACES " \t\r\n"

namespace gstlrn
{
class AMatrix;

#ifndef SWIG
Id _getColumnRank();
Id _getColumnName();
Id _getColumnSize(Id localSize = 0);
Id _getMaxNCols();
Id _getMaxNRows();
Id _getNBatch();
Id _getDecimalNumber();
double _getThresh();

String _toStrRowColumn(Id icase, Id value, Id flagAdd);
String _toStrTrailer(Id ncols, Id nrows, Id ncols_util, Id nrows_util);
String _toStrColumnHeader(const VectorString& colnames, Id colfrom, Id colto, Id colSize = 0);
String _toStrRowHeader(const VectorString& rownames, Id iy, Id rowSize = 0);
#endif

GSTLEARN_EXPORT String toUpper(const std::string_view string);
GSTLEARN_EXPORT String toLower(const std::string_view string);
#ifndef SWIG
GSTLEARN_EXPORT void toUpper(String& string);
GSTLEARN_EXPORT void toLower(String& string);
#endif // SWIG

GSTLEARN_EXPORT void skipBOM(std::ifstream& ins);
GSTLEARN_EXPORT bool matchKeyword(const String& string1,
                                  const String& string2,
                                  bool caseSensitive = true);
GSTLEARN_EXPORT bool matchRegexp(const String& string1,
                                 const String& string2,
                                 bool caseSensitive = true);
GSTLEARN_EXPORT Id getRankInList(const VectorString& list,
                                 const String& match,
                                 bool caseSensitive = true);
GSTLEARN_EXPORT Id decodeInString(const String& symbol,
                                  const String& node,
                                  Id* facies,
                                  bool caseSensitive = true);
GSTLEARN_EXPORT Id decodeInList(const VectorString& symbols,
                                const String& node,
                                Id* rank,
                                Id* facies,
                                bool caseSensitive = true);
GSTLEARN_EXPORT void correctNamesForDuplicates(VectorString& newList,
                                               const VectorString& oldList = VectorString());
GSTLEARN_EXPORT void correctNewNameForDuplicates(VectorString& list, Id target);

GSTLEARN_EXPORT String incrementStringVersion(const String& string,
                                              Id rank             = 1,
                                              const String& delim = ".");
GSTLEARN_EXPORT String concatenateString(const String& string,
                                         double value,
                                         const String& delim = "-");
GSTLEARN_EXPORT String concatenateStrings(const String& delimt  = ".",
                                          const String& string1 = "",
                                          const String& string2 = "",
                                          const String& string3 = "",
                                          const String& string4 = "");

GSTLEARN_EXPORT VectorString generateMultipleNames(const String& radix,
                                                   Id number,
                                                   const String& delim = "-");
GSTLEARN_EXPORT VectorString expandList(const VectorString& list,
                                        const String& match,
                                        bool onlyOne = false);
GSTLEARN_EXPORT VectorString expandList(const VectorString& list,
                                        const VectorString& matches);
GSTLEARN_EXPORT Id getMaxStringSize(const VectorString& list);
GSTLEARN_EXPORT VectorString separateKeywords(const String& code);

GSTLEARN_EXPORT String trimRight(const String& s, const String& t = SPACES);
GSTLEARN_EXPORT String trimLeft(const String& s, const String& t = SPACES);
GSTLEARN_EXPORT String trim(const String& s, const String& t = SPACES);
GSTLEARN_EXPORT String erase(const String& s, const String& t = SPACES);

GSTLEARN_EXPORT VectorInt decodeGridSorting(const String& string,
                                            const VectorInt& nx,
                                            bool verbose = false);

/**
 * @brief Convert the contents of any argument (double, Id, String) into a String
 *
 * @tparam T Can be double, Id or String
 * @param v Identified argument
 * @param justification -1 for Left justified; 0 for center; 1 for right justification
 * @param localSize Dimension provided for the formatted output string
 * @return String Returned string
 */
template<typename T>
String toStr(const T& v, Id justification = 1, Id localSize = 0);
String toStr(const char* v, Id justification = 1, Id localSize = 0);
String toStrVectorVec(const String& title, constvect tab, bool flagOverride = true);
String toStrTitle(Id level, const char* format, ...);
String toStrInterval(double zmin, double zmax);
VectorString toStrVectorDouble(const VectorDouble& values, Id justification = 1);

/**
 * @brief Converting the contents of a String into double or integer
 *
 * @tparam T Must be a String
 * @param v Input string
 * @param dec Number of decimals (only used for Double conversion)
 * @return T Returned vlaue (double or Id)
 */
template<typename T>
T fromStr(const String& v, char dec = '.');

/**
 * @brief Ask for a value interactively
 *
 * @tparam T Input string
 * @param v Input string
 * @param defval Default value
 * @param authTest Is an NA answer authorized
 * @return T Valid for Double, Id or Bool
 */
template<typename T>
T askInteractive(const String& v, T defval, bool authTest = false);

/**
 * Print the contents of a VectorDouble in a Matrix Form
 * @param title        Title of the printout
 * @param colnames     Names of the columns (optional)
 * @param rownames     Names of the rows (optional)
 * @param bycol        True if values as sorted by column; false otherwise
 * @param nrows        Number of rows
 * @param ncols        Number of columns
 * @param tab          VectorDouble containing the values
 * @param flagOverride True to override printout limitations
 * @param flagSkipZero when True, skip the zero values (represented by a '.' as for sparse matrix)
 */
template<typename T>
inline String toStrMatrix(const String& title,
                          const VectorString& colnames,
                          const VectorString& rownames,
                          bool bycol,
                          Id nrows,
                          Id ncols,
                          const VectorNumT<T>& tab,
                          bool flagOverride = false,
                          bool flagSkipZero = false)
{
  // Vérification du type à la compilation
  static_assert(std::is_same<T, Id>::value ||
                  std::is_same<T, double>::value,
                "toStrMatrix: T must be int or double");
  std::stringstream sstr;
  if (tab.empty() || ncols <= 0 || nrows <= 0) return sstr.str();

  /* Initializations */

  Id ncutil = ncols;
  Id nrutil = nrows;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();
  Id npass       = static_cast<Id>(ceil(static_cast<double>(ncutil) / static_cast<double>(_getNBatch())));
  bool multi_row = nrutil > 1 || npass > 1;

  Id colSize = 0;
  if (colnames.empty())
    colSize = _getColumnSize();
  else
    colSize = MAX(MIN(_getColumnName(), getMaxStringSize(colnames) + 1), _getColumnSize());
  Id rowSize = 0;
  if (rownames.empty())
    rowSize = _getColumnSize();
  else
    rowSize = MAX(getMaxStringSize(rownames) + 1, _getColumnSize());

  /* Print the title (optional) */

  if (!title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  // Loop on the batches

  for (Id ipass = 0; ipass < npass; ipass++)
  {
    Id jdeb = ipass * _getNBatch();
    Id jfin = MIN(jdeb + _getNBatch(), ncutil);

    /* Print the names of the columns and the column numbers */

    if (multi_row) sstr << _toStrColumnHeader(colnames, jdeb, jfin, colSize);

    /* Loop on the rows */

    for (Id iy = 0; iy < nrutil; iy++)
    {
      if (multi_row) sstr << _toStrRowHeader(rownames, iy, rowSize);

      /* Loop on the columns */
      for (Id ix = jdeb; ix < jfin; ix++)
      {
        Id iad = (bycol) ? iy + nrows * ix : ix + ncols * iy;
        if (flagSkipZero && ABS(tab[iad]) < EPSILON20)
          sstr << toStr(".", 1, colSize);
        else
          sstr << toStr(tab[iad], 1, colSize);
      }
      sstr << std::endl;
    }
  }

  /* Print the trailer */

  sstr << _toStrTrailer(ncols, nrows, ncutil, nrutil);
  return sstr.str();
}
GSTLEARN_EXPORT String toMatrix(const String& title,
                                const AMatrix& mat,
                                bool flagOverride = false,
                                bool flagSkipZero = false);

template<typename T>
inline String toStrVector(const String& title,
                          const VectorT<T>& tab,
                          bool flagOverride = true)
{
  static_assert(std::is_same<T, Id>::value ||
                  std::is_same<T, double>::value ||
                  std::is_same<T, String>::value,
                "toStrVector: T must be Id, double, or String");

  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  Id ncols  = static_cast<Id>(tab.size());
  Id ncutil = ncols;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagOverride) ncutil = _getMaxNCols();
  bool multi_row = ncutil > _getNBatch();

  /* Print the title (optional) */

  if (!title.empty())
  {
    sstr << title;
    if (multi_row) sstr << std::endl;
  }

  Id lec = 0;
  if (multi_row) sstr << _toStrColumnHeader(VectorString(), 0, _getNBatch());

  for (Id i = 0; i < ncutil; i += _getNBatch())
  {
    if (multi_row) sstr << _toStrRowHeader(VectorString(), i);

    for (Id j = 0; j < _getNBatch(); j++)
    {
      if (lec >= ncutil) continue;
      sstr << toStr(tab[lec]);
      lec++;
    }
    sstr << std::endl;
  }

  // Print the trailer
  sstr << _toStrTrailer(ncols, 0, ncutil, 0);

  return sstr.str();
}
template<typename T>
inline String toStrVector(const String& title,
                          const VectorT<VectorNumT<T>>& tab,
                          bool flagOverride = true)
{
  static_assert(std::is_same<T, Id>::value ||
                  std::is_same<T, double>::value,
                "toStrVector: T must be Id or double");
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  if (!title.empty())
    sstr << title << std::endl;

  Id nrows  = static_cast<Id>(tab.size());
  Id nrutil = nrows;
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagOverride) nrutil = _getMaxNRows();

  for (Id i = 0; i < nrutil; i++)
    sstr << toStrVector(String(), tab[i], flagOverride);

  // Print the trailer
  sstr << _toStrTrailer(0, nrows, 0, nrutil);

  return sstr.str();
}

// Functions calling modern tools and regrouping potential warnings issued by Windows
GSTLEARN_EXPORT char* gslStrcpy(char* dst, Id n, const char* src);
GSTLEARN_EXPORT void gslStrcpy(String& dst, const String& src);
GSTLEARN_EXPORT void gslStrcpy(String& dst, const char* src);
GSTLEARN_EXPORT char* gslStrcat(char* dst, Id n, const char* src);
GSTLEARN_EXPORT void gslStrcat(String& dst, const char* src);
GSTLEARN_EXPORT void gslStrcat(String& dst, const String& src);
GSTLEARN_EXPORT Id gslSPrintf(char* dst, Id n, const char* fmt, ...);
GSTLEARN_EXPORT Id gslSPrintf(String& dst, const char* fmt, ...);
GSTLEARN_EXPORT Id gslSPrintfCat(String& dst, const char* fmt, ...);
GSTLEARN_EXPORT char* gslStrtok(char* str, const char* delim);
GSTLEARN_EXPORT Id gslScanf(const char* fmt, ...);
GSTLEARN_EXPORT Id gslSScanf(const char* str, const char* fmt, ...);

} // namespace gstlrn