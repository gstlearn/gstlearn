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

#include "Basic/StringHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

// TODO : add Namespace
#define SPACES " \t\r\n"

namespace gstlrn
{
class AMatrix;

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
 * @brief Ask for a value interactively
 *
 * @param v Input string
 * @param defval Default value
 * @param authTest If an NA answer is authorized
 */
// Template générique (déclaré mais pas défini)
template<typename T>
T question(const String& v, T defval = getNA<T>(), bool authTest = false);

// Spécialisation pour double
template<>
inline double question<double>(const String& v, double defval, bool authTest)
{
  return _askDouble(v, defval, authTest);
}

// Spécialisation pour Id
template<>
inline Id question<Id>(const String& v, Id defval, bool authTest)
{
  return _askInt(v, defval, authTest);
}

// Spécialisation pour bool
template<>
inline bool question<bool>(const String& v, bool defval, bool authTest)
{
  return _askBool(v, defval, authTest);
}

/**
 * @brief Convert the contents of any argument (double, Id, String) into a String
 *
 * @param v Identified argument
 * @param justification -1 for Left justified; 0 for center; 1 for right justification
 * @param localSize Dimension provided for the formatted output string
 * @param roundZero true to round very small double values to zero
 * @param nColumns Number of columns for the formatted output string
 * @param flagScientific true to use scientific notation for double values
 * @return String Returned string
 */
GSTLEARN_EXPORT String toStr(double v,
                             Id justification    = 1,
                             Id localSize        = 0,
                             bool roundZero      = true,
                             Id nColumns         = 1,
                             bool flagScientific = false);
GSTLEARN_EXPORT String toStr(Id v,
                             Id justification = 1,
                             Id localSize     = 0,
                             bool roundZero   = true,
                             Id nColumns      = 1);
GSTLEARN_EXPORT String toStr(const String& v,
                             Id justification = 1,
                             Id localSize     = 0,
                             bool roundZero   = true,
                             Id nColumns      = 1);
GSTLEARN_EXPORT String toStrFormat(const char* format, ...);

GSTLEARN_EXPORT String toStrTitle(Id level, const char* format, ...);
GSTLEARN_EXPORT String toStrInterval(double zmin, double zmax);
GSTLEARN_EXPORT String toStrVectorVec(const String& title,
                                      constvect tab,
                                      bool flagIgnoreMaxNCols = false,
                                      bool newLineAfterTitle  = true);

/**
 * @brief Converting the contents of a String into double or integer
 *
 * @param s Input string
 * @param dec Number of decimals (only used for Double conversion)
 */
// Template générique (non défini)
template<typename T>
T fromStr(const String& s, char dec = '.');

// Spécialisation pour double
template<>
inline double fromStr<double>(const String& s, char dec)
{
  return _convertToDouble(s, dec);
}

// Spécialisation pour Id
template<>
inline Id fromStr<Id>(const String& s, char dec)
{
  return _convertToId(s, dec);
}

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
String toStrMatrix(const String& title,
                   const VectorString& colnames,
                   const VectorString& rownames,
                   bool bycol,
                   Id nrows,
                   Id ncols,
                   const VectorNumT<T>& tab,
                   bool flagOverride = false,
                   bool flagSkipZero = false)
{
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
    if (multi_row) sstr << _toStrColumnHeaders(colnames, jdeb, jfin, colSize);

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

GSTLEARN_EXPORT String toStrMatrix(const String& title,
                                   const AMatrix& mat,
                                   bool flagOverride = false,
                                   bool flagSkipZero = false);
GSTLEARN_EXPORT String toStrTrimMat(const VectorDouble& tl,
                                    Id neq,
                                    Id mode             = 1,
                                    const String& title = String());

/**
 * Printout a vector in a formatted manner
 * @param title Title of the printout
 * @param tab   Vector of values to be printed
 * @param flagIgnoreMaxNCols true to ignore the maximum number of columns
 * @param newLineAfterTitle true to put a new line after the title (if any)
 * @return The string (terminated with a newline)
 */
template<typename T>
String toStrVector(const String& title,
                   const VectorT<T>& tab,
                   bool flagIgnoreMaxNCols = false,
                   bool newLineAfterTitle  = false)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  Id ncols  = static_cast<Id>(tab.size());
  Id ncutil = ncols;
  if (_getMaxNCols() > 0 && ncutil > _getMaxNCols() && !flagIgnoreMaxNCols) ncutil = _getMaxNCols();
  bool multi_row = ncols > _getNBatch();

  /* Print the title (optional) */

  if (!title.empty())
  {
    sstr << title;
    if (newLineAfterTitle) sstr << std::endl;
  }

  Id lec = 0;
  if (multi_row) sstr << _toStrColumnHeaders(VectorString(), 0, _getNBatch());

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

/**
 * Printout a vector in a formatted manner
 * @param title Title of the printout
 * @param tab  VectorVectorX to be printed
 * @param flagIgnoreMaxNRows true to ignore the maximum number of rows
 * @param newLineAfterTitle true to put a new line after the title (if any)
 * @return The string (terminated with a newline)
 */
template<typename T>
String toStrVector(const String& title,
                   const VectorNumT<VectorNumT<T>>& tab,
                   bool flagIgnoreMaxNRows = false,
                   bool newLineAfterTitle  = true)
{
  std::stringstream sstr;
  if (tab.empty()) return sstr.str();

  if (!title.empty())
  {
    sstr << title;
    if (newLineAfterTitle) sstr << std::endl;
  }

  Id nrows  = static_cast<Id>(tab.size());
  Id nrutil = nrows;
  if (_getMaxNRows() > 0 && nrutil > _getMaxNRows() && !flagIgnoreMaxNRows) nrutil = _getMaxNRows();

  for (Id i = 0; i < nrutil; i++)
    sstr << toStrVector(String(), tab[i], flagIgnoreMaxNRows);

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