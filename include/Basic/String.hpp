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
 * @brief Ask for a value interactively
 *
 * @param v Input string
 * @param defval Default value
 * @param authTest If an NA answer is authorized
 */
GSTLEARN_EXPORT double questionDouble(const String& v, double defval = TEST, bool authTest = false);
GSTLEARN_EXPORT Id questionId(const String& v, Id defval = ITEST, bool authTest = false);
GSTLEARN_EXPORT bool questionBool(const String& v, bool defval = false, bool authTest = false);

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
 * @param string Input string
 * @param dec Number of decimals (only used for Double conversion)
 */
GSTLEARN_EXPORT double fromStrToDouble(const String& string, char dec = '.');
GSTLEARN_EXPORT Id fromStrToId(const String& string, char dec = '.');

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
                   bool flagSkipZero = false);

extern template String toStrMatrix<double>(const String& title,
                                           const VectorString& colnames,
                                           const VectorString& rownames,
                                           bool bycol,
                                           Id nrows,
                                           Id ncols,
                                           const VectorNumT<double>& tab,
                                           bool flagOverride = false,
                                           bool flagSkipZero = false);
extern template String toStrMatrix<long long>(const String& title,
                                              const VectorString& colnames,
                                              const VectorString& rownames,
                                              bool bycol,
                                              Id nrows,
                                              Id ncols,
                                              const VectorNumT<long long>& tab,
                                              bool flagOverride = false,
                                              bool flagSkipZero = false);

GSTLEARN_EXPORT String toMatrix(const String& title,
                                const AMatrix& mat,
                                bool flagOverride = false,
                                bool flagSkipZero = false);
GSTLEARN_EXPORT String toStrTrimMat(const String& title,
                                    Id mode,
                                    Id neq,
                                    const VectorDouble& tl);

/**
 * Printout a template vector in a formatted manner
 * @param title Title of the printout
 * @param tab   Template (VectorAny) to be printed
 * @param flagIgnoreMaxNCols true to ignore the maximum number of columns
 * @param newLineAfterTitle true to put a new line after the title (if any)
 * @return The string (terminated with a newline)
 */
template<typename T>
String toStrVector(const String& title,
                   const VectorT<T>& tab,
                   bool flagIgnoreMaxNCols = false,
                   bool newLineAfterTitle  = false);

extern template String toStrVector<double>(const String& title,
                                           const VectorT<double>& tab,
                                           bool flagIgnoreMaxNCols,
                                           bool newLineAfterTitle);
extern template String toStrVector<long long>(const String& title,
                                              const VectorT<long long>& tab,
                                              bool flagIgnoreMaxNCols,
                                              bool newLineAfterTitle);
extern template String toStrVector<String>(const String& title,
                                           const VectorT<String>& tab,
                                           bool flagIgnoreMaxNCols,
                                           bool newLineAfterTitle);

/**
 * Printout a template vector in a formatted manner
 * @param title Title of the printout
 * @param tab   Template (VectorVectorDouble or VectorVectorInt) to be printed
 * @param flagIgnoreMaxNRows true to ignore the maximum number of rows
 * @param newLineAfterTitle true to put a new line after the title (if any)
 * @return The string (terminated with a newline)
 */
template<typename T>
String toStrVector(const String& title,
                   const VectorNumT<VectorNumT<T>>& tab,
                   bool flagIgnoreMaxNRows = false,
                   bool newLineAfterTitle  = true);

extern template String toStrVector<double>(const String& title,
                                           const VectorNumT<VectorNumT<double>>& tab,
                                           bool flagIgnoreMaxNRows,
                                           bool newLineAfterTitle);
extern template String toStrVector<long long>(const String& title,
                                              const VectorNumT<VectorNumT<long long>>& tab,
                                              bool flagIgnoreMaxNRows,
                                              bool newLineAfterTitle);

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