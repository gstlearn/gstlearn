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
#include "math.h"
#include "Basic/String.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include <algorithm>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <regex>

/**
 * Protect the matching pattern against Crash which happens when the string
 * contains "*" without any preceding character
 * @param match Initial matching pattern
 * @return The std::regex item used for further comparisons
 */
std::regex _protectRegexp(const String& match)
{
  String str = match;
  std::size_t found = str.find('*');
  if (found != std::string::npos)
  {
    if (found == 0 || str[found-1] != '.') str.insert(found,".");
  }
  std::regex regexpr(str);

  return regexpr;
}

String toUpper(const String& string)
{
  String str = string;
  toUpper(str);
  return (str);
}

String toLower(const String& string)
{
  String str = string;
  toLower(str);
  return (str);
}

void toUpper(String& string)
{
  std::for_each(string.begin(), string.end(), [](char & c){
      c = ::toupper(c);
  });
}

void toLower(String& string)
{
  std::for_each(string.begin(), string.end(), [](char & c){
      c = ::tolower(c);
  });
}

enum charTypeT{ other, alpha, digit};

charTypeT _charType(char c){
    if(isdigit(c))return digit;
    if(isalpha(c))return alpha;
    return other;
}

String incrementStringVersion(const String& string,
                              int rank,
                              const String& delim)
{
  std::stringstream ss;
  ss << string << delim << rank;
  return ss.str();
}

String concatenateStrings(const String& delim,
                          const String& string1,
                          const String& string2,
                          const String& string3,
                          const String& string4)
{
  bool started = false;
  std::stringstream ss;
  if (! string1.empty())
  {
    started = true;
    ss << string1;
  }
  if (! string2.empty())
  {
    if (started) ss << delim;
    ss << string2;
    started = true;
  }
  if (! string3.empty())
  {
    if (started) ss << delim;
    ss << string3;
    started = true;
  }
  if (! string4.empty())
  {
    if (started) ss << delim;
    ss << string4;
    started = true;
  }
  return ss.str();
}

VectorString generateMultipleNames(const String& radix,
                                   int number)
{
  VectorString list;

  for (int i = 0; i < number; i++)
  {
    list.push_back(incrementStringVersion(radix, i+1));
  }
  return list;
}

/**
 * Check that the name 'rank' is not conflicting with any previous name.
 * If it does, increment its name by a version number.
 * @param list
 * @return Number of items whose names are modified
 */
int correctNamesForDuplicates(VectorString& list)
{
  int numberRenamed = 0;
  int number = static_cast<int> (list.size());
  for (int i = 1; i < number; i++)
  {
    // Check that a similar name does not appear among the previous names in list

    int found = -1;
    for (int j = 0; j < i && found < 0; j++)
    {
      if (list[i].compare(list[j]) == 0) found = j;
    }
    if (found < 0) continue;

    // We have found a similar name

    while (list[found].compare(list[i]) == 0)
      list[i] = incrementStringVersion(list[i]);
    numberRenamed++;
  }
  return numberRenamed;
}

void correctNewNameForDuplicates(VectorString& list, int rank)
{
   int number = static_cast<int> (list.size());
   int found = 1;
   while (found > 0)
   {
     found = 0;
     for (int i = 0; i < number; i++)
     {
       if (i == rank) continue;
       if (list[rank].compare(list[i]) == 0) found++;
     }
     if (found <= 0) break;;

     // We have found a similar name

     list[rank] = incrementStringVersion(list[rank]);
   }
}

/**
 * Return the rank of the (first) item in 'list' which matches 'match'
 * @param list  List of keywords used for search
 * @param match Searched pattern
 * @param caseSensitive Case Sensitive flag
 * @return The index of the matching item or -1
 */
int getRankInList(const VectorString& list,
                  const String& match,
                  bool caseSensitive)
{
  for (int i = 0; i< (int) list.size(); i++)
  {
    if (matchRegexp(list[i], match, caseSensitive)) return i;
  }
  return -1;
}

/**
 * Decode the input string 'node' as a one-character keyword followed by an integer rank
 * The keyword must match the symbol.
 * The integer rank is returned as 'facies'
 * @param symbol
 * @param node
 * @param facies
 * @param caseSensitive
 * @return Error returned code
 */
int decodeInString(const String& symbol,
                   const String& node,
                   int *facies,
                   bool caseSensitive)
{
  String locnode = node;
  String locsymb = symbol;
  if (!caseSensitive)
  {
    toUpper(locnode);
    toUpper(locsymb);
  }

  if (locnode.compare(0, 1, locsymb) != 0) return 1;

  // Decode the remaining part: it should be an integer

  for (char &c : locnode)
  {
    if (!isdigit(c)) c = ' ';
  }

  *facies = 0;
  std::stringstream ss(locnode);
  ss >> *facies;
  return 0;
}

/**
 * Decode the input string 'node' as a one-character keyword followed by an integer rank
 * The keyword must match one of the symbols: its rank is 'rank'
 * The integer rank is returned as 'facies'
 * @param symbols
 * @param node
 * @param rank
 * @param facies
 * @param caseSensitive
 * @return Error returned code
 */
int decodeInList(const VectorString& symbols,
                 const String& node,
                 int *rank,
                 int *facies,
                 bool caseSensitive)
{
  String local = node;

  for (int i = 0; i < (int) symbols.size(); i++)
  {
    if (decodeInString(symbols[i], node, facies, caseSensitive)) continue;
    *rank = i;
    return 0;
  }
  return -1;
}

/**
 * Check if two keywords match up to "Regular Expression" on the second string and their case
 * @param string1 First keyword
 * @param string2 Second keyword (through Regular Expression interpretation)
 * @param caseSensitive true if the case of both strings should be considered
 * Otherwise, both strings are converted into upper case before comparison
 * @return true if both keywords are identical; false otherwise
 */
bool matchRegexp(const String& string1, const String& string2, bool caseSensitive)
{
  String local1 = string1;
  String local2 = string2;
  if (! caseSensitive)
  {
    toUpper(local1);
    toUpper(local2);
  }
  std::regex regexpr = _protectRegexp(local2);
  if (std::regex_match(local1, regexpr)) return true;
  return false;
}

/**
 * Check if two keywords are similar, up to their case
 * @param string1 First keyword
 * @param string2 Second keyword
 * @param caseSensitive true if the case of both strings should be considered
 * Otherwise, both strings are converted into upper case before comparison
 * @return true if both keywords are identical; false otherwise
 */
bool matchKeyword(const String& string1, const String& string2, bool caseSensitive)
{
  String local1 = string1;
  String local2 = string2;
  if (! caseSensitive)
  {
    toUpper(local1);
    toUpper(local2);
  }
  return local1 == local2;
}

/**
 * Returns the list of matching names
 * @param list   List of eligible names
 * @param match  Name to be expanded
 * @param onlyOne True if the expanded list may only contain a single name
 * @return The list of matching possibilities
 * @note The returned list is empty if no match has been found or if
 *       several matches have been found but onlyOne flag is True (message issued).
 * @remark Example:  expandList(["x", "y", "xcol"], "x.*") -> ("x", "xcol")
 */
VectorString expandList(const VectorString& list,
                        const String& match,
                        bool onlyOne)
{
  std::regex regexpr = _protectRegexp(match);

  VectorString sublist;
  for (int i = 0; i< (int) list.size(); i++)
  {
    if (std::regex_match(list[i], regexpr))
      sublist.push_back(list[i]);
  }

  int number = static_cast<int> (sublist.size());
  if (onlyOne && number != 1)
  {
    if (number > 1)
    {
      messerr("The name (%s) has been expanded to several matching possibilities",
              match.c_str());
      for (int i = 0; i < (int) sublist.size(); i++)
        messerr("- %s",sublist[i].c_str());
    }
    else
      messerr("The name (%s) does not have any matching possibility",
              match.c_str());
    messerr("A single match is requested");
    sublist.clear();
  }
  return sublist;
}

VectorString expandList(const VectorString& list,
                        const VectorString& matches)
{
  VectorString sublist;

  // Loop on the eligible names
  for (int i = 0; i< (int) list.size(); i++)
  {
    // Loop on the target names
    bool found = false;
    for (int j = 0; j < (int) matches.size() && ! found; j++)
    {
      std::regex regexpr = _protectRegexp(matches[j]);
      if (std::regex_match(list[i], regexpr)) found = true;
    }
    if (found) sublist.push_back(list[i]);
  }
  return sublist;
}

/**
 * Returns the maximum string size of a list of strings
 * @param list List of strings
 * @return The maximum number of characters
 */
int getMaxStringSize(const VectorString& list)
{
  int size = 0;
  if (list.empty()) return size;
  for (int i = 0; i < (int) list.size(); i++)
  {
    int local = static_cast<int> (list[i].length());
    if (local > size) size = local;
  }
  return size;
}

VectorString separateKeywords(const String& code)
{
  VectorString result;
  String oString = "";
  charTypeT st = other;
  for (auto c : code)
  {
    if ((st == alpha && _charType(c) == digit) || (st == digit && _charType(c) == alpha))
    {
      if (oString.size() > 0) result.push_back(oString);
      oString = "";
    }
    oString.push_back(c);
    st = _charType(c);
  }
  if (oString.size() > 0) result.push_back(oString);
  return result;
}

/**
 * Decode an integer from a string. Returns ITEST if impossible
 * @param code String to be decoded
 * @return The integer value or ITEST (in case of failure)
 */
int toInt(const String& code)
{
  std::istringstream iss(code);
  int number;
  iss >> number;
  if (iss.fail())
    return ITEST;
  else
    return number;
}

String intToString(int value)
{
  std::stringstream sstr;
  sstr << value;
  return sstr.str();
}

String realToString(double value)
{
  std::stringstream sstr;
  sstr << value;
  return sstr.str();
}

String suppressTrailingBlanks(String value)
{
  value = std::regex_replace(value, std::regex(" +$"), "");
  return value;
}

String suppressLeadingBlanks(String value)
{
  value = std::regex_replace(value, std::regex("^ +"), "");
  return value;
}

String suppressAnyBlanks(String value)
{
  value = std::regex_replace(value, std::regex("^ +| +$|( ) +"), "$1");
  return value;
}

/**
 * Ask interactively for the value of one integer
 * @param text Text of the question
 * @param defval Default value (or IFFFF)
 */
int askInt(const String& text, int defval)
{
  bool nodef = ! IFFFF(defval);
  int answer = defval;
  std::cin.exceptions(std::istream::failbit|
                      std::istream::badbit);

  try
  {
    String str;
    do {
      if (nodef)
        std::cout << text << " (Default = " << defval << ") : ";
      else
        std::cout << text << " : ";
      std::getline(std::cin, str);
    }
    while (! nodef && str.empty());
    std::stringstream sstr(str);
    sstr >> answer; // TODO verifier que le resultat et bien un entier
  }
  catch(std::istream::failure e)
  {
    std::cerr << "Problem when reading integer" << std::endl;
  }
  return answer;
}

/**
 * Ask interactively for the value of one Real (Double)
 * @param text Text of the question
 * @param defval Default value (or IFFFF)
 */
double askDouble(const String& text, double defval)
{
  bool nodef = ! FFFF(defval);
  double answer = defval;
  std::cin.exceptions(std::istream::failbit|
                      std::istream::badbit);

  try
  {
    String str;
    do {
      if (nodef)
        std::cout << text << " (Default = " << defval << ") : ";
      else
        std::cout << text << " : ";
      std::getline(std::cin, str);
    }
    while (! nodef && str.empty());
    std::stringstream sstr(str);
    sstr >> answer; // TODO: verifier que le resultat est bien un reel, sinon???
  }
  catch(std::istream::failure e)
  {
    std::cerr << "Problem when reading integer" << std::endl;
  }
  return answer;
}
