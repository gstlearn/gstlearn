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
#pragma once

#include "geoslib_define.h"
#include "Basic/AException.hpp"

#include <memory>
#include <string>

// TODO : add Namespace
#define SPACES " \t\r\n"

void skipBOM(std::ifstream &in);

String toUpper(const String& string);
String toLower(const String& string);

void toUpper(String& string);
void toLower(String& string);

bool matchKeyword(const String& string1,
                  const String& string2,
                  bool caseSensitive = true);
bool matchRegexp(const String& string1,
                 const String& string2,
                 bool caseSensitive = true);
int getRankInList(const VectorString& list,
                  const String& string,
                  bool caseSensitive = true);
int decodeInString(const String& symbol,
                   const String& node,
                   int *facies,
                   bool caseSensitive = true);
int decodeInList(const VectorString& symbols,
                 const String& node,
                 int *rank,
                 int *facies,
                 bool caseSenstive = true);
int  correctNamesForDuplicates(VectorString& list);
void correctNewNameForDuplicates(VectorString& list, int rank);

String incrementStringVersion(const String& string,
                              int rank = 1,
                              const String& delim = ".");
String concatenateStrings(const String& delimt = ".",
                          const String& string1 = String(),
                          const String& string2 = String(),
                          const String& string3 = String(),
                          const String& string4 = String());

VectorString generateMultipleNames(const String& radix, int number);
VectorString expandList(const VectorString& list,
                        const String& match,
                        bool onlyOne = false);
VectorString expandList(const VectorString& list,
                        const VectorString& matches);
int getMaxStringSize(const VectorString& list);
VectorString separateKeywords(const String& code);

// TODO : Use template functions
int toInt(const String& v);
double toDouble(const String& v, char dec = '.');
String toString(int value);
String toString(double value);
int    askInt(   const String& text, int    defval = ITEST, bool authTest = false);
double askDouble(const String& text, double defval = TEST,  bool authTest = false);
int    askBool(  const String& text, bool   defval = false);

String trimRight (const String& s, const String& t = SPACES);
String trimLeft(  const String& s, const String& t = SPACES);
String trim(      const String& s, const String& t = SPACES);
String erase(     const String& s, const String& t = SPACES);

char* gslStrcpy(char* dst, const char* src);
char* gslStrcat(char* dst, const char* src);
int   gslSPrintf(char* dst, const char* fmt, ...);

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

