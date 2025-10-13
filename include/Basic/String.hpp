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
GSTLEARN_EXPORT void skipBOM(std::ifstream& ins);

GSTLEARN_EXPORT String toUpper(const std::string_view string);
GSTLEARN_EXPORT String toLower(const std::string_view string);

#ifndef SWIG
GSTLEARN_EXPORT void toUpper(String& string);
GSTLEARN_EXPORT void toLower(String& string);
#endif // SWIG

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
GSTLEARN_EXPORT void correctNamesForDuplicates(VectorString& list);
GSTLEARN_EXPORT void correctNewNameForDuplicates(VectorString& list, Id rank);

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

// TODO : Use template functions
GSTLEARN_EXPORT Id toInteger(const String& v);
GSTLEARN_EXPORT double toDouble(const String& v, char dec = '.');
GSTLEARN_EXPORT String toString(Id value);
GSTLEARN_EXPORT String toString(double value);
GSTLEARN_EXPORT Id askInt(const String& text,
                          Id defval     = ITEST,
                          bool authTest = false);
GSTLEARN_EXPORT double askDouble(const String& text,
                                 double defval = TEST,
                                 bool authTest = false);
GSTLEARN_EXPORT Id askBool(const String& text, bool defval = false);

GSTLEARN_EXPORT String trimRight(const String& s, const String& t = SPACES);
GSTLEARN_EXPORT String trimLeft(const String& s, const String& t = SPACES);
GSTLEARN_EXPORT String trim(const String& s, const String& t = SPACES);
GSTLEARN_EXPORT String erase(const String& s, const String& t = SPACES);

GSTLEARN_EXPORT VectorInt decodeGridSorting(const String& string,
                                            const VectorInt& nx,
                                            bool verbose = false);

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