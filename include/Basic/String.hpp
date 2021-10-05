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

#include "Basic/AStringable.hpp"

// TODO : add Namespace

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
int toInt(const String& code);
String intToString(int value);
String realToString(double value);
// TODO : add const ref
String suppressTrailingBlanks(String value);
String suppressLeadingBlanks(String value);
String suppressAnyBlanks(String value);


