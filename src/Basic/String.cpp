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
#include "Basic/String.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/VectorNumT.hpp"
#include "Matrix/AMatrix.hpp"

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
// Functions for returning parameters from OptCst environment
Id _getColumnRank()
{
  return static_cast<Id>(OptCst::query(ECst::NTRANK));
}
Id _getColumnName()
{
  return static_cast<Id>(OptCst::query(ECst::NTNAME));
}
Id _getColumnSize(Id localSize)
{
  return (localSize > 0) ? localSize : static_cast<Id>(OptCst::query(ECst::NTCAR));
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

// Functions for encoding a single value
std::stringstream _formatOneColumn(Id justification, Id localSize = 0)
{
  std::stringstream sstr;
  auto size = static_cast<I32>(_getColumnSize(localSize));
  auto prec = static_cast<I32>(_getDecimalNumber());
  sstr << std::fixed << std::setw(size) << std::setprecision(prec);
  if (justification < 0)
    sstr << std::left;
  else
    sstr << std::right;
  return sstr;
}

String _formatOneString(const String& string, Id justification, Id localSize)
{
  std::stringstream sstr = _formatOneColumn(justification, localSize);
  Id size                = static_cast<Id>(string.size());
  Id truncSize           = _getColumnSize(localSize);
  if (size > truncSize)
  {
    // String must be truncated

    String strloc = string;
    strloc.erase(0, size - truncSize);
    strloc.replace(0, 2, " *");
    sstr << strloc;
  }
  else
  {
    sstr << string;
  }
  return sstr.str();
}

String _formatOneDouble(double value, Id justification, Id localSize)
{
  std::stringstream sstr = _formatOneColumn(justification, localSize);
  if (FFFF(value))
    sstr << "N/A";
  else
  {
    // Prevent -0.00 : https://stackoverflow.com/a/12536500/3952924
    value = (ABS(value) < _getThresh()) ? 0. : value;
    sstr << value;
  }

  return sstr.str();
}

String _formatOneId(Id value, Id justification, Id localSize = 0)
{
  std::stringstream sstr = _formatOneColumn(justification, localSize);
  if (IFFFF(value))
    sstr << "N/A";
  else
    sstr << value;

  return sstr.str();
}

/**
 * Decode an double from a string. Returns TEST if impossible
 * @param v String to be decoded
 * @param dec Decimal separator character
 * @return The double value or TEST (in case of failure)
 */
template<typename T>
class dec_separator: public std::numpunct<T>
{
public:
  dec_separator(char dec = ',')
    : _dec(dec)
  {
  }

private:
  typename std::numpunct<T>::char_type do_decimal_point() const override
  {
    return _dec;
  }
  char _dec;
};

/**
 * Decode an integer from a string. Returns ITEST when impossible
 * @param v String to be decoded
 * @param dec Symbol used for decimal point (unused here)
 * @return The integer value or ITEST (in case of failure)
 */
Id _convertToInteger(const String& v, char dec)
{
  DECLARE_UNUSED(dec);
  std::istringstream iss(v);
  Id number;
  iss >> number;
  if (iss.fail()) return ITEST;
  return number;
}

/**
 * Decode a double value from a string. Returns TEST when impossible
 * @param v String to be decoded
 * @param dec Symbol used for decimal point
 * @return The double value or TEST (in case of failure)
 */
double _convertToDouble(const String& v, char dec)
{
  std::istringstream iss(v);
  double number;
  iss.imbue(std::locale(iss.getloc(), new dec_separator<char>(dec)));
  iss >> number;
  if (iss.fail()) return TEST;
  return number;
}

/**
 * Ask interactively for the value of one integer
 * @param v Text of the question
 * @param defval Default value (or IFFFF)
 * @param authTest True if TEST value is authorized (TEST)
 */
Id _askInt(const String& v, Id defval, bool authTest)
{
  bool hasDefault = !IFFFF(defval) || authTest;
  Id answer       = defval;
  std::cin.exceptions(std::istream::failbit | std::istream::badbit);

  try
  {
    while (true)
    {
      // Display the question
      if (hasDefault)
      {
        if (IFFFF(defval))
          std::cout << v << " (Default = TEST) : ";
        else
          std::cout << v << " (Default = " << defval << ") : ";
      }
      else
        std::cout << v << " : ";

      // Read the answer
      String str;
      std::getline(std::cin, str);

      // Check for empty line: set to default value
      if (str.empty() && hasDefault)
      {
        answer = defval;
        break;
      }

      // Check the TEST answer

      if (authTest && str == "TEST")
      {
        answer = ITEST;
        break;
      }

      // Try casting in integer
      std::stringstream ss(str);
      if (ss >> answer) break;

      std::cout << "The answer is not a valid integer!" << std::endl;
    }
  }
  catch (std::istream::failure& e)
  {
    std::cerr << "Problem when reading integer:" << e.what() << std::endl;
  }
  return answer;
}

/**
 * Ask interactively for the value of one Real (Double value)
 * @param v Text of the question
 * @param defval Default value (or IFFFF)
 * @param authTest True if a TEST answer is authorized (TEST)
 */
double _askDouble(const String& v, double defval, bool authTest)
{
  bool hasDefault = !FFFF(defval) || authTest;
  double answer   = defval;
  std::cin.exceptions(std::istream::failbit | std::istream::badbit);

  try
  {
    while (true)
    {
      // Display the question
      if (hasDefault)
      {
        if (FFFF(defval))
          std::cout << v << " (Default = TEST) : ";
        else
          std::cout << v << " (Default = " << defval << ") : ";
      }
      else
        std::cout << v << " : ";

      // Read the answer
      String str;
      std::getline(std::cin, str);

      // Check for empty line: set to default value
      if (str.empty() && hasDefault)
      {
        answer = defval;
        break;
      }

      // Catch the TEST answer
      if (authTest && str == "TEST")
      {
        answer = TEST;
        break;
      }

      // Try casting in integer
      std::stringstream ss(str);
      if (ss >> answer) break;

      std::cout << "The answer is not a valid double!" << std::endl;
    }
  }
  catch (std::istream::failure& e)
  {
    std::cerr << "Problem when reading double:" << e.what() << std::endl;
  }
  return answer;
}

/**
 * Ask interactively for the value of one boolean
 * @param v Text of the question
 * @param defval Default value
 * @param authTest True if a TEST answer is authorized (TEST)
 */
bool _askBool(const String& v, bool defval, bool authTest)
{
  DECLARE_UNUSED(authTest);
  bool hasDefault = !IFFFF(defval);
  bool answer     = defval;
  std::cin.exceptions(std::istream::failbit | std::istream::badbit);

  try
  {
    while (true)
    {
      // Display the question
      if (hasDefault)
      {
        String defstr;
        if (defval)
          defstr = "Y";
        else
          defstr = "N";
        std::cout << v << " (Default = " << defstr << ") : ";
      }
      else
        std::cout << v << " : ";

      // Read the answer
      String str;
      std::getline(std::cin, str);

      // Check for empty line: set to default value
      if (str.empty() && hasDefault)
      {
        answer = defval;
        break;
      }

      // Try checking authorized answer
      if (str == "Y")
      {
        answer = true;
        break;
      }
      if (str == "N")
      {
        answer = false;
        break;
      }

      std::cout << "The answer is not a valid bool!" << std::endl;
    }
  }
  catch (std::istream::failure& e)
  {
    std::cerr << "Problem when reading bool:" << e.what() << std::endl;
  }
  return answer;
}
template<>
String toStr(const double& v, Id justification, Id localSize)
{
  std::stringstream sstr;
  sstr << _formatOneDouble(v, justification, localSize);
  return sstr.str();
}
template<>
String toStr(const Id& v, Id justification, Id localSize)
{
  std::stringstream sstr;
  sstr << _formatOneId(v, justification, localSize);
  return sstr.str();
}
template<>
String toStr(const String& v, Id justification, Id localSize)
{
  std::stringstream sstr;
  sstr << _formatOneString(v, justification, localSize);
  return sstr.str();
}

String toStr(const char* v, Id justification, Id localSize)
{
  return toStr(String(v), justification, localSize);
}

/**
 * Print a message and underlines it with various formats
 * @param level  Level of the title
 * @param format Output format
 * @param ...    Additional arguments
 */
String toStrTitle(Id level, const char* format, ...)
{
  std::stringstream sstr;
  Id STRING_MAX = 1000;
  char STRING[1000];
  va_list ap;

  sstr << std::endl;
  va_start(ap, format);
  (void)vsnprintf(STRING, sizeof(STRING), format, ap);
  va_end(ap);
  sstr << STRING << std::endl;

  /* Underline the string */

  Id size = static_cast<Id>(strlen(STRING));
  (void)gslStrcpy(STRING, STRING_MAX, "");
  for (Id i = 0; i < size; i++)
  {
    switch (level)
    {
      case 0:
        (void)gslStrcat(STRING, STRING_MAX, "=");
        break;

      case 1:
        (void)gslStrcat(STRING, STRING_MAX, "-");
        break;

      case 2:
        (void)gslStrcat(STRING, STRING_MAX, ".");
        break;
    }
  }
  sstr << STRING << std::endl;

  return sstr.str();
}

String toStrInterval(double zmin, double zmax)
{
  std::stringstream sstr;

  sstr << "[";
  if (FFFF(zmin))
    sstr << "N/A";
  else
    sstr << zmin;
  message(" ; ");
  if (FFFF(zmax))
    sstr << "N/A";
  else
    sstr << zmax;
  sstr << "]" << std::endl;

  return sstr.str();
}

/**
 * Printout a vector in a formatted manner
 * @param title Title of the printout (or empty string)
 * @param tab   Vector (real values) to be printed
 * @param flagOverride true to override printout limitations
 * @return The string (terminated with a newline)
 */
String toStrVectorVec(const String& title, constvect tab, bool flagOverride)
{
  // Adapter le span vers un VectorT<double> ou équivalent
  VectorT<double> tmp(tab.begin(), tab.end());

  return toStrVector(title, tmp, flagOverride);
}

VectorString toStrVectorDouble(const VectorDouble& values, Id justification)
{
  VectorString strings;
  for (Id i = 0; i < static_cast<Id>(values.size()); i++)
    strings.push_back(toStr(values[i], justification));
  return strings;
}

// Functions for beautifying a suite of values
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

String _toStrRowHeader(const VectorString& rownames, Id iy, Id rowSize)
{
  std::stringstream sstr;
  if (!rownames.empty())
    sstr << toStr(rownames[iy], -1, rowSize);
  else
    sstr << _toStrRowColumn(CASE_ROW, iy, false);
  return sstr.str();
}
String _toStrColumnHeader(const VectorString& colnames, Id colfrom, Id colto, Id colSize)
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

/**
 * Protect the matching pattern against Crash which happens when the string
 * contains "*" without any preceding character
 * @param match Initial matching pattern
 * @return The std::regex item used for further comparisons
 */
std::regex _protectRegexp(const String& match)
{
  String str        = match;
  std::size_t found = str.find('*');
  if (found != String::npos)
  {
    if (found == 0 || str[found - 1] != '.') str.insert(found, ".");
  }
  std::regex regexpr(str);

  return regexpr;
}

String toUpper(const std::string_view string)
{
  String str {string};
  toUpper(str);
  return (str);
}

String toLower(const std::string_view string)
{
  String str {string};
  toLower(str);
  return (str);
}

void toUpper(String& string)
{
  std::for_each(string.begin(), string.end(), [](char& c)
                { c = static_cast<char>(::toupper(c)); });
}

void toLower(String& string)
{
  std::for_each(string.begin(), string.end(), [](char& c)
                { c = static_cast<char>(::tolower(c)); });
}

enum charTypeT
{
  other,
  alpha,
  digit
};

charTypeT _charType(char c)
{
  if (isdigit(c)) return digit;
  if (isalpha(c)) return alpha;
  return other;
}

String incrementStringVersion(const String& string,
                              Id rank,
                              const String& delim)
{
  std::stringstream ss;
  ss << string << delim << rank;
  return ss.str();
}

String concatenateString(const String& string,
                         double value,
                         const String& delim)
{
  std::stringstream ss;
  ss << trim(string) << delim << value;
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
  if (!string1.empty())
  {
    started = true;
    ss << trim(string1);
  }
  if (!string2.empty())
  {
    if (started) ss << delim;
    ss << trim(string2);
    started = true;
  }
  if (!string3.empty())
  {
    if (started) ss << delim;
    ss << trim(string3);
    started = true;
  }
  if (!string4.empty())
  {
    if (started) ss << delim;
    ss << trim(string4);
  }
  return ss.str();
}

VectorString generateMultipleNames(const String& radix, Id number, const String& delim)
{
  VectorString list;

  for (Id i = 0; i < number; i++)
  {
    list.push_back(incrementStringVersion(radix, i + 1, delim));
  }
  return list;
}

/**
 * Check that the names in 'list' are not conflicting with any previous name.
 * If it does, increment its name by a version number.
 * @param newList List of new names (suggested in Input and possibly corrected)
 * @param oldList List of already accepted names
 */
void correctNamesForDuplicates(VectorString& newList,
                               const VectorString& oldList)
{
  Id nnew = static_cast<Id>(newList.size());
  Id nold = static_cast<Id>(oldList.size());

  VectorString catList = oldList;                                // start with contents of 'a'
  catList.insert(catList.end(), newList.begin(), newList.end()); // append contents of 'b'

  for (Id i = 0; i < nnew; i++)
  {
    // Check that a similar name does not appear among the previous names in list
    String nameref = newList[i];
    if (nameref.empty()) continue;

    Id rank = 0;
  label_try:
    rank++;
    Id found = -1;
    for (Id j = 0; j < i + nold && found < 0; j++)
    {
      if (newList[i] == catList[j]) found = j;
    }
    if (found < 0) continue;

    // We have found a similar name. Modify it as long as it matches an already existing name

    newList[i] = incrementStringVersion(nameref, rank);
    goto label_try;
  }
}

void correctNewNameForDuplicates(VectorString& list, Id itarget)
{
  Id number      = static_cast<Id>(list.size());
  Id found       = 1;
  String nameref = list[itarget];
  while (found > 0)
  {
    Id rank = 0;

  label_try:
    rank++;
    found = -1;
    for (Id i = 0; i < number && found < 0; i++)
    {
      if (i == itarget) continue;
      if (list[itarget] == list[i]) found++;
    }
    if (found < 0) break;

    // We have found a similar name

    list[itarget] = incrementStringVersion(nameref, rank);
    goto label_try;
  }
}

/**
 * Return the rank of the (first) item in 'list' which matches 'match'
 * @param list  List of keywords used for search
 * @param match Searched pattern
 * @param caseSensitive Case Sensitive flag
 * @return The index of the matching item or -1
 */
Id getRankInList(const VectorString& list,
                 const String& match,
                 bool caseSensitive)
{
  for (Id i = 0; i < static_cast<Id>(list.size()); i++)
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
Id decodeInString(const String& symbol,
                  const String& node,
                  Id* facies,
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

  for (char& c: locnode)
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
Id decodeInList(const VectorString& symbols,
                const String& node,
                Id* rank,
                Id* facies,
                bool caseSensitive)
{
  for (Id i = 0; i < static_cast<Id>(symbols.size()); i++)
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
bool matchRegexp(const String& string1,
                 const String& string2,
                 bool caseSensitive)
{
  String local1 = string1;
  String local2 = string2;
  if (!caseSensitive)
  {
    toUpper(local1);
    toUpper(local2);
  }
  std::regex regexpr = _protectRegexp(local2);
  return std::regex_match(local1, regexpr);
}

/**
 * Check if two keywords are similar, up to their case
 * @param string1 First keyword
 * @param string2 Second keyword
 * @param caseSensitive true if the case of both strings should be considered
 * Otherwise, both strings are converted into upper case before comparison
 * @return true if both keywords are identical; false otherwise
 */
bool matchKeyword(const String& string1,
                  const String& string2,
                  bool caseSensitive)
{
  String local1 = string1;
  String local2 = string2;
  if (!caseSensitive)
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
  for (Id i = 0; i < static_cast<Id>(list.size()); i++)
  {
    const String& toto = list[i];
    if (std::regex_match(toto, regexpr)) sublist.push_back(toto);
  }

  Id number = static_cast<Id>(sublist.size());
  if (onlyOne && number != 1)
  {
    if (number > 1)
    {
      messerr(
        "The name (%s) has been expanded to several matching possibilities",
        match.c_str());
      for (Id i = 0; i < static_cast<Id>(sublist.size()); i++)
        messerr("- %s", sublist[i].c_str());
    }
    else
      messerr("The name (%s) does not have any matching possibility",
              match.c_str());
    messerr("A single match is requested");
    sublist.clear();
  }
  return sublist;
}

VectorString expandList(const VectorString& list, const VectorString& matches)
{
  VectorString sublist;

  // Loop on the patterns to be matched
  for (Id i = 0; i < static_cast<Id>(matches.size()); i++)
  {
    // Loop for eligible names
    for (Id j = 0; j < static_cast<Id>(list.size()); j++)
    {
      std::regex regexpr = _protectRegexp(matches[i]);
      if (std::regex_match(list[j], regexpr))
      {
        // A match is found

        // Check that the name has not been recorded yet in sublist
        bool already = false;
        if (std::find(sublist.begin(), sublist.end(), list[j]) != sublist.end())
        {
          already = true;
        }

        // If the name is not already registered: add it
        if (!already) sublist.push_back(list[j]);
      }
    }
  }
  return sublist;
}

/**
 * Returns the maximum string size of a list of strings
 * @param list List of strings
 * @return The maximum number of characters
 */
Id getMaxStringSize(const VectorString& list)
{
  Id size = 0;
  if (list.empty()) return size;
  for (Id i = 0; i < static_cast<Id>(list.size()); i++)
  {
    Id local = static_cast<Id>(list[i].length());
    if (local > size) size = local;
  }
  return size;
}

/**
 * Separate keywords in the input string
 * The keywords are identified when switching between alpha, digit and other
 * @param code String to be split
 * @return
 */
VectorString separateKeywords(const String& code)
{
  VectorString result;
  String oString;
  charTypeT st = other;
  for (auto c: code)
  {
    if ((st == alpha && _charType(c) != alpha) || (st == digit && _charType(c) != digit) || (st == other && _charType(c) != other))
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

String trimRight(const String& s, const String& t)
{
  String d(s);
  String::size_type i(d.find_last_not_of(t));
  if (i == String::npos) return "";
  return d.erase(d.find_last_not_of(t) + 1);
}

String trimLeft(const String& s, const String& t)
{
  String d(s);
  return d.erase(0, s.find_first_not_of(t));
}

String trim(const String& s, const String& t)
{
  const String& d(s);
  return trimLeft(trimRight(d, t), t);
}

String erase(const String& s, const String& t)
{
  String d(s);
  for (size_t i = 0; i < t.size(); i++)
    d.erase(std::remove(d.begin(), d.end(), t[i]), d.end());
  return d;
}

char* gslStrcpy(char* dst, Id n, const char* src)
{
  return strncpy(dst, src, n);
}

void gslStrcpy(String& dst, const char* src)
{
  if (!src)
  {
    dst.clear();
    return;
  }

  size_t len = std::strlen(src);

  // ajuster la taille du string
  dst.resize(len);

  // copier le contenu
  std::memcpy(dst.data(), src, len);
}

void gslStrcpy(String& dst, const String& src)
{
  size_t len = src.size();

  // redimensionner dst pour accueillir exactement src
  dst.resize(len);

  // copier le contenu
  std::memcpy(dst.data(), src.data(), len);
}

char* gslStrcat(char* dst, Id n, const char* src)
{
  return strncat(dst, src, n);
}

void gslStrcat(String& dst, const char* src)
{
  if (!src) return; // sécurité

  size_t old_len = dst.size();
  size_t add_len = std::strlen(src);

  // on redimensionne pour accueillir l'ajout
  dst.resize(old_len + add_len);

  // on copie le nouveau contenu à la fin
  std::memcpy(&dst[old_len], src, add_len);
}

void gslStrcat(String& dst, const String& src)
{
  size_t old_len = dst.size();
  size_t add_len = src.size();

  // redimensionner dst pour accueillir l'ajout
  dst.resize(old_len + add_len);

  // copier directement les caractères de src
  std::memcpy(&dst[old_len], src.data(), add_len);
}

Id gslSPrintf(char* dst, Id n, const char* fmt, ...)
{
  DECLARE_UNUSED(n);
  va_list ap;
  va_start(ap, fmt);
  Id nout = vsnprintf(dst, n, fmt, ap);
  va_end(ap);
  return nout;
}

Id gslSPrintfCat(String& dst, const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);

  // calculer la taille du texte formaté
  va_list ap_copy;
  va_copy(ap_copy, ap);
  int size = std::vsnprintf(nullptr, 0, fmt, ap_copy);
  va_end(ap_copy);

  if (size < 0)
  {
    va_end(ap);
    return -1; // erreur
  }

  // mémoriser l'ancienne longueur
  size_t old_len = dst.size();

  // redimensionner dst pour contenir l'ancien contenu + nouveau texte + '\0'
  dst.resize(old_len + size + 1);

  // écrire le texte formaté à la fin
  std::vsnprintf(&dst[old_len], size + 1, fmt, ap);
  va_end(ap);

  // redimensionner pour enlever le '\0' final
  dst.resize(old_len + size);

  return size;
}

Id gslSPrintf(String& dst, const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);

  // On calcule la taille nécessaire
  va_list ap_copy;
  va_copy(ap_copy, ap);
  int size = std::vsnprintf(nullptr, 0, fmt, ap_copy);
  va_end(ap_copy);

  if (size < 0)
  {
    va_end(ap);
    return -1; // erreur
  }

  // Redimensionner pour contenir le texte + '\0'
  dst.resize(size + 1);

  // Ecrire la chaîne formatée
  std::vsnprintf(dst.data(), dst.size(), fmt, ap);
  va_end(ap);

  // Enlever le '\0' final pour que dst.size() == longueur réelle
  dst.resize(size);

  return size;
}

Id gslScanf(const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  Id n = vscanf(fmt, ap);
  va_end(ap);
  return n;
}

Id gslSScanf(const char* str, const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  Id n = vsscanf(str, fmt, ap);
  va_end(ap);
  return n;
}

char* gslStrtok(char* str, const char* delim)
{
  return strtok(str, delim);
}

/****************************************************************************/
/*!
 **  Decode the grid sorting order
 **
 ** \return Array describing the order
 **
 ** \param[in]  string  Name of the sorting string
 ** \param[in]  nx      Array giving the number of cells per direction
 ** \param[in]  verbose Verbose flag
 **
 ** \remarks The value of order[i] gives the dimension of the space along
 ** \remarks which sorting takes place at rank "i". This value is positive
 ** \remarks for increasing order and negative for decreasing order.
 ** \remarks 'order' values start from 1 (for first space dimension)
 ** \remarks Example: "+x1-x3+x2"
 **
 *****************************************************************************/
VectorInt decodeGridSorting(const String& string,
                            const VectorInt& nx,
                            bool verbose)
{
  Id ndim = static_cast<Id>(nx.size());
  VectorInt order(ndim, 0);
  VectorInt ranks(ndim, 0);

  // Loop on the character string

  Id idim   = 0;
  Id ind    = 0;
  Id length = static_cast<Id>(string.size());

  while (ind < length)
  {
    Id orient = 0;
    if (string[ind] == '-' && (string[ind + 1] == 'x' && isdigit(string[ind + 2])))
      orient = -1;
    else if (string[ind] == '+' && (string[ind + 1] == 'x' && isdigit(string[ind + 2])))
      orient = 1;
    else
      orient = 0;

    if (orient != 0)
    {

      // The string '+x*' or '-x*' has been encountered

      ind += 2;
      Id num = string[ind] - '0';
      if (idim >= ndim)
      {
        messerr("'order' contains more terms (%d) than the space dimension (%d)", idim + 1, ndim);
        return VectorInt();
      }
      order[idim] = orient * num;
      if (num > ndim)
      {
        messerr("'order' refers to 'x%d' while space dimension is %d", num, ndim);
        return VectorInt();
      }
      ranks[num - 1] = orient;
      idim++;
    }
    else
    {
      ind++;
    }
  }

  // Check that all indices (within space dimension) have been specified

  for (Id i = 0; i < ndim; i++)
  {
    if (ranks[i] != 0) continue;
    messerr("'x%d' is not mentioned in the input string", i + 1);
    return VectorInt();
  }

  // Optional printout

  if (verbose)
  {
    message("Decoding the sorting rule (%s) with nx = (", string.c_str());
    for (Id i = 0; i < ndim; i++)
      message(" %d", nx[i]);
    message(" )\n");
    for (Id i = 0; i < ndim; i++)
    {
      Id a_order = ABS(order[i]);
      message("%d - Dimension=%d - N%d=%d", i + 1, a_order, a_order,
              nx[a_order - 1]);
      if (order[i] > 0)
        message(" - Increasing\n");
      else
        message(" - Decreasing\n");
    }
  }
  return (order);
}

template<>
double fromStr(const String& v, char dec)
{
  return _convertToDouble(v, dec);
}
template<>
Id fromStr(const String& v, char dec)
{
  return _convertToInteger(v, dec);
}

template<>
double askInteractive(const String& v, double defval, bool authTest)
{
  return _askDouble(v, defval, authTest);
}
template<>
Id askInteractive(const String& v, Id defval, bool authTest)
{
  return _askInt(v, defval, authTest);
}
template<>
bool askInteractive(const String& v, bool defval, bool authTest)
{
  return _askBool(v, defval, authTest);
}

/**
 * @overload
 * Print the contents of a Matrix in a Matrix Form
 * @param title        Title of the printout
 * @param mat          Contents of a AMatrix
 * @param flagOverride true to override printout limitations
 * @param flagSkipZero when true, skip the zero values (represented by a '.' as for sparse matrix)
 *                     always true for sparse matrix
 */
String toMatrix(const String& title,
                const AMatrix& mat,
                bool flagOverride,
                bool flagSkipZero)
{
  flagSkipZero = mat.isSparse();
  return toStrMatrix(title, VectorString(), VectorString(), true,
                     mat.getNRows(), mat.getNCols(), mat.getValues(),
                     flagOverride, flagSkipZero);
}

} // namespace gstlrn