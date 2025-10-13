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
#include "OutputFormat/FileLAS.hpp"
#include "Basic/String.hpp"
#include "Db/Db.hpp"
#include "OutputFormat/AOF.hpp"

#include <cctype>
#include <cstring>

#if defined(_WIN32) || defined(_WIN64)
namespace
{
// STL replacement for strcasestr - far from optimal (several copies)
// C++23 introduces std::string::contains
// TODO: should be moved elsewhere...
std::string tolower(const std::string& s)
{
  std::string result;
  result.reserve(s.size());
  std::transform(s.begin(), s.end(), std::back_inserter(result),
                 [](unsigned char c)
                 { return std::tolower(c); });
  return result;
};
const char* strcasestr(const char* s, const char* ss)
{
  const auto pos = tolower(s).find(tolower(ss));
  if (pos == std::string::npos) return (const char*)nullptr;
  return s + pos;
};
} // namespace
#endif

namespace gstlrn
{
FileLAS::FileLAS(const char* filename, const Db* db)
  : AOF(filename, db)
  , _xwell(0.)
  , _ywell(0.)
  , _cwell(0.)
{
}

FileLAS::FileLAS(const FileLAS& r)
  : AOF(r)
  , _xwell(r._xwell)
  , _ywell(r._ywell)
  , _cwell(r._cwell)
{
}

FileLAS& FileLAS::operator=(const FileLAS& r)
{
  if (this != &r)
  {
    AOF::operator=(r);
    _xwell = r._xwell;
    _ywell = r._ywell;
    _cwell = r._cwell;
  }
  return *this;
}

FileLAS::~FileLAS()
{
}

Db* FileLAS::readFromFile()
{
  Db* db = nullptr;
  Id string_max = 1000;
  char string[1000], *lcur, sep_blank[2], sep_point[2], *token;
  double value;
  static I32 s_length = 1000;
  VectorString names  = {"X", "Y", "CODE"};

  /* Open the file */

  if (_fileReadOpen()) return db;

  // Initializations */

  double test  = TEST;
  sep_blank[0] = ' ';
  sep_blank[1] = '\0';
  sep_point[0] = '.';
  sep_point[1] = '\0';
  Id nvar      = static_cast<Id>(names.size());

  // Loop on the lines

  Id numline = 0;
  (void)gslStrcpy(string, string_max, "");

  // Decode the header
  if (_readFind(s_length, "~Well", &numline, string)) return db;
  while (1)
  {
    if (_readNext(s_length, 1, &numline, string)) return db;

    // Looking for the TEST value (keyword "NULL")
    if (strstr(string, "NULL")) gslSScanf(&string[8], "%lf", &test);

    // Looking for the next delimitor
    if (strstr(string, "~")) break;
  }

  // Decode the variable names
  if (_readFind(s_length, "~Curve", &numline, string)) return db;
  while (1)
  {
    if (_readNext(s_length, 1, &numline, string)) return db;

    // Skipping the comments
    if (strstr(string, "#")) continue;

    // Looking for the variable header
    if (strstr(string, "~")) break;

    // Reading the variable name
    token = gslStrtok(string, sep_point);
    if (token == NULL) break;

    // Add the variable to the list

    string[strlen(token)] = '\0';
    names.push_back(string);
    nvar++;
  }

  /* Decoding the array of data */

  if (_readFind(s_length, "~A", &numline, string)) return db;

  VectorDouble tab;
  Id nech = 0;
  while (1)
  {
    if (_readNext(s_length, 1, &numline, string)) break;

    // Add the special variables

    tab.push_back(_xwell);
    tab.push_back(_ywell);
    tab.push_back(_cwell);
    Id nvarlu = 3;

    // Add other variables

    lcur = string;
    while (nvarlu < nvar)
    {
      token = gslStrtok(lcur, sep_blank);
      if (token == nullptr) break;
      lcur = NULL;
      if (gslSScanf(token, "%lf", &value) == EOF) break;
      if (value == test) value = TEST;
      tab.push_back(value);
      nvarlu++;
    }

    if (nvarlu < nvar)
    {
      if (fgets(string, s_length, _file) == nullptr) break;
    }
    nech++;
  }

  db = new Db();
  db->resetFromSamples(nech, ELoadBy::SAMPLE, tab, names);

  // Close the file

  _fileClose();

  return db;
}

/****************************************************************************/
/*!
 **  Read the next lines from a LAS file until the search file is found
 **
 ** \return  0 if the string is found; 1 for end-of-file
 **
 ** \param[in]      s_length   Length of the string array
 ** \param[in]      target     Target line to be found
 ** \param[in,out]  numline    Rank of the line
 ** \param[out]     string     New line
 **
 *****************************************************************************/
Id FileLAS::_readFind(I32 s_length,
                      const char* target,
                      Id* numline,
                      char* string)
{
  /* Check the current line */
  if (strcasestr(string, target)) return (0);

  while (1)
  {
    if (_readNext(s_length, 1, numline, string)) return (1);
    if (strcasestr(string, target)) return (0);
  }
  return (1);
}

/****************************************************************************/
/*!
 **   Read the next line from a LAS file
 **
 ** \return  Error return code
 **
 ** \param[in]      s_length   Length of the string array
 ** \param[in]      flag_up    Convert to uppercase
 ** \param[in,out]  numline    Rank of the line
 ** \param[out]     string     New line
 **
 ** TODO : Replace char* by String
 **
 *****************************************************************************/
Id FileLAS::_readNext(I32 s_length, Id flag_up, Id* numline, char* string)
{
  Id size;

  (*numline)++;
  if (fgets(string, s_length, _file) == nullptr) return (1);
  size = static_cast<Id>(strlen(string));

  // Suppress the trailing newline
  if (string[size - 1] == '\n') string[size - 1] = '\0';
  size = static_cast<Id>(strlen(string));

  // Replace any illegal character from the string by a blank
  for (Id i = 0; i < size; i++)
  {
    if (string[i] < ' ' || string[i] > '~') string[i] = ' ';
  }

  // Convert to uppercase (optional)
  if (flag_up) _stringToUppercase(string);

  return (0);
}

void FileLAS::_stringToUppercase(char* string)
{
  Id n = static_cast<Id>(strlen(string));
  for (Id i = 0; i < n; i++)
    if (string[i] >= 'a' && string[i] <= 'z')
      string[i] = ('A' + string[i] - 'a');
}
} // namespace gstlrn
