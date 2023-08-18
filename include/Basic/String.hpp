/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include "Basic/AException.hpp"
#include "Basic/VectorT.hpp"
#include "Basic/VectorNumT.hpp"

#include <memory>

// TODO : add Namespace
#define SPACES " \t\r\n"

GSTLEARN_EXPORT void skipBOM(std::ifstream &ins);

GSTLEARN_EXPORT String toUpper(const String &string);
GSTLEARN_EXPORT String toLower(const String &string);

GSTLEARN_EXPORT void toUpper(String &string);
GSTLEARN_EXPORT void toLower(String &string);

GSTLEARN_EXPORT bool matchKeyword(const String &string1,
                                  const String &string2,
                                  bool caseSensitive = true);
GSTLEARN_EXPORT bool matchRegexp(const String &string1,
                                 const String &string2,
                                 bool caseSensitive = true);
GSTLEARN_EXPORT int getRankInList(const VectorString &list,
                                  const String &string,
                                  bool caseSensitive = true);
GSTLEARN_EXPORT int decodeInString(const String &symbol,
                                   const String &node,
                                   int *facies,
                                   bool caseSensitive = true);
GSTLEARN_EXPORT int decodeInList(const VectorString &symbols,
                                 const String &node,
                                 int *rank,
                                 int *facies,
                                 bool caseSenstive = true);
GSTLEARN_EXPORT void correctNamesForDuplicates(VectorString &list);
GSTLEARN_EXPORT void correctNewNameForDuplicates(VectorString &list, int rank);

GSTLEARN_EXPORT String incrementStringVersion(const String &string,
                                              int rank = 1,
                                              const String &delim = ".");
GSTLEARN_EXPORT String concatenateString(const String &string,
                                         double value,
                                         const String &delim = "-");
GSTLEARN_EXPORT String concatenateStrings(const String &delimt = ".",
                                          const String &string1 = String(),
                                          const String &string2 = String(),
                                          const String &string3 = String(),
                                          const String &string4 = String());

GSTLEARN_EXPORT VectorString generateMultipleNames(const String &radix,
                                                   int number,
                                                   const String& delim = "-");
GSTLEARN_EXPORT VectorString expandList(const VectorString &list,
                                        const String &match,
                                        bool onlyOne = false);
GSTLEARN_EXPORT VectorString expandList(const VectorString &list,
                                        const VectorString &matches);
GSTLEARN_EXPORT int getMaxStringSize(const VectorString &list);
GSTLEARN_EXPORT VectorString separateKeywords(const String &code);

// TODO : Use template functions
GSTLEARN_EXPORT int toInteger(const String &v);
GSTLEARN_EXPORT double toDouble(const String &v, char dec = '.');
GSTLEARN_EXPORT String toString(int value);
GSTLEARN_EXPORT String toString(double value);
GSTLEARN_EXPORT int askInt(const String &text,
                           int defval = ITEST,
                           bool authTest = false);
GSTLEARN_EXPORT double askDouble(const String &text,
                                 double defval = TEST,
                                 bool authTest = false);
GSTLEARN_EXPORT int askBool(const String &text, bool defval = false);

GSTLEARN_EXPORT String trimRight(const String &s, const String &t = SPACES);
GSTLEARN_EXPORT String trimLeft(const String &s, const String &t = SPACES);
GSTLEARN_EXPORT String trim(const String &s, const String &t = SPACES);
GSTLEARN_EXPORT String erase(const String &s, const String &t = SPACES);

GSTLEARN_EXPORT VectorInt decodeGridSorting(const String& name,
                                            const VectorInt& nx,
                                            bool verbose = false);

GSTLEARN_EXPORT char* gslStrcpy(char *dst, const char *src);
GSTLEARN_EXPORT char* gslStrcat(char *dst, const char *src);
GSTLEARN_EXPORT int gslSPrintf(char *dst, const char *fmt, ...);
GSTLEARN_EXPORT char* gslStrtok(char *str, const char *delim);
GSTLEARN_EXPORT char* gslStrncpy(char *dest, const char *src, size_t n);
GSTLEARN_EXPORT int gslScanf(const char *format, ...);
GSTLEARN_EXPORT int gslSScanf(const char *str, const char *format, ...);
GSTLEARN_EXPORT int gslFScanf(FILE *stream, const char *format, ...);

// Adapted from:
// - https://stackoverflow.com/a/26310318
// - https://stackoverflow.com/a/26221725

/**
 * Secured version of sprintf (using char *)
 *
 * @param dst Destination string
 * @param fmt Formatted string
 * @param args Additional arguments
 *
 * @remark: dst must be pre allocated with appropriate size
 */
/*
 template<typename... Args>
 int gslSPrintf(char* dst, const char* fmt, Args... args)
 {
 size_t size_s = std::snprintf(nullptr, 0, fmt, args...);
 if (size_s == 0) { my_throw("Error during formatting."); }
 snprintf(dst, size_s + 1, fmt, args...);
 return static_cast<int>(size_s + 1);
 }
 */
/**
 * Secured version of sprintf (using String)
 *
 * @param dst Destination string
 * @param fmt Formatted string
 * @param args Additional arguments
 *
 */
/*
 template<typename ... Args>
 int gslSPrintf(String& dst, String fmt, Args... args )
 {
 int size_s = std::snprintf( nullptr, 0, fmt.c_str(), args ... ) + 1; // Extra space for '\0'
 if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
 size_t size = static_cast<size_t>( size_s );
 char* buf = new char(size_s); // make_unique not yet available in c++11
 std::snprintf( buf, size, fmt.c_str(), args ... );
 dst = std::string( buf, buf + size - 1 ); // We don't want the '\0' inside
 delete buf;
 return static_cast<int>(size_s + 1);
 }
 */

