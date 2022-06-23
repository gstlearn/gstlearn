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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Basic/AStringable.hpp"
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"

#include <iostream>
#include <stdarg.h>
#include <stdio.h>
#include <fstream>

class GSTLEARN_EXPORT ASerializable
{
public:
  ASerializable();
  ASerializable(const ASerializable& r);
  ASerializable& operator=(const ASerializable& r);
  virtual ~ASerializable();

  bool deserialize(std::istream& is, bool verbose = false);
  bool serialize(std::ostream& os,bool verbose = false) const;
  bool dumpToNF(const String& neutralFilename, bool verbose = false) const;

  static String buildFileName(const String& filename, bool ensureDirExist = false);

  static String getHomeDirectory(const String& sub = "");
  static String getWorkingDirectory();
  static String getTestData(const String& subdir, const String& filename);
  static String getFileIdentity(const String& filename);
  static void setContainerName(bool useDefault,
                               const String& containerName = String(),
                               bool verbose = false);
  static void unsetContainerName();
  static void setPrefixName(const String& prefixName);
  static void unsetPrefixName();
  static const String& getContainerName();
  static const String& getPrefixName();

  // TODO : Directory manipulation class
  static bool createDirectory(const String& dir);
  static String getExecDirectory();
  static String getDirectory(const String& path);

protected:
  virtual bool _deserialize(std::istream& s, bool verbose = false) = 0;
  virtual bool _serialize(std::ostream& s,bool verbose = false) const = 0;
  virtual String _getNFName() const = 0;

  static bool _fileOpenWrite(const String& filename,
                              const String& filetype,
                              std::ofstream& os,
                              bool verbose = false);
  static bool _fileOpenRead(const String& filename,
                             const String& filetype,
                             std::ifstream& is,
                             bool verbose = false);

  static bool _commentWrite(std::ostream& os,
                             const String& comment);
  template <typename T>
  static bool _recordWrite(std::ostream& os,
                            const String& title,
                            const T& val);
  template <typename T>
  static bool _recordWriteVec(std::ostream& os,
                               const String& title,
                               const std::vector<T>& vec);

  template <typename T>
  static bool _recordRead(std::istream& is,
                           const String& title,
                           T& val);
  template <typename T>
  static bool _recordReadVec(std::istream& is,
                              const String& title,
                              std::vector<T>& vec,
                              int nvalues);

  static bool _onlyBlanks(char *string);

  static int  _tableRead(std::istream& is, int ntab, double *tab);
  static bool _tableWrite(std::ostream& os,
                           const String& string,
                           int ntab,
                           const VectorDouble& tab);

private:
  static String myContainerName;
  static String myPrefixName;
};

template <typename T>
bool ASerializable::_recordWrite(std::ostream& os,
                                 const String& title,
                                 const T& val)
{
  if (os.good())
  {
    if (isNA<T>(val))
    {
      if (title.empty())
         os << STRING_NA << " ";
       else
         os << STRING_NA << " # " << title << std::endl;
    }
    else
    {
      if (title.empty())
        os << val << " ";
      else
        os << val << " # " << title << std::endl;
    }
  }
  return os.good();
}

template <typename T>
bool ASerializable::_recordWriteVec(std::ostream& os,
                                    const String& title,
                                    const std::vector<T>& vec)
{
  if (os.good())
  {
    if (!title.empty())
      os << "# " << title << std::endl;
    for (auto val: vec)
    {
      if (isNA<T>(val))
        os << STRING_NA << " ";
      else
        os << val << " ";
    }
    os << std::endl;
  }
  return os.good();
}

template <typename T>
bool ASerializable::_recordRead(std::istream& is, const String& title, T& val)
{
  val = T();
  if (is.good())
  {
    String word;
    // Skip comment or empty lines
    while (is.good())
    {
      word.clear();
      is >> word;
      if (!is.good() && !is.eof())
      {
        messerr("Error while reading %s", title.c_str());
        return false;
      }
      word = trim(word);
      if (!word.empty())
      {
        if (word == STRING_NA) break;   // We found NA
        else if (word[0] != '#') break; // We found something
        else std::getline(is, word);    // We found comment, eat all the line
      }
    }

    if (word == STRING_NA)
    {
      // Get NA value
      val = getNAValue<T>();
    }
    else
		{
      // Decode the line
      std::stringstream sstr(word);
      sstr >> val;
      if (!sstr.good() && !sstr.eof())
      {
        messerr("Error while reading %s", title.c_str());
        val = T();
        return false;
      }
    }
  }
  return true;
}

template <typename T>
bool ASerializable::_recordReadVec(std::istream& is,
                                   const String& title,
                                   std::vector<T>& vec,
                                   int nvalues)
{
  vec.clear();
  if (is.good())
  {
    String line;
    // Skip comment or empty lines
    while (is.good())
    {
      std::getline(is, line);
      if (!is.good() && !is.eof())
      {
        messerr("Error while reading %s", title.c_str());
        return false;
      }
      line = trim(line);
      if (!line.empty() && line[0] != '#')
        break; // We found something
    }
    // Decode the line
    std::stringstream sstr(line);
    while (sstr.good())
    {
      String word;
      sstr >> word;
      if (!sstr.good() && !sstr.eof())
      {
        messerr("Error while reading %s", title.c_str());
        vec.clear();
        return false;
      }
      word = trim(word);
      if (word.empty()) continue;
      T val;
      if (word == STRING_NA)
      {
        // Get NA value
        val = getNAValue<T>();
      }
      else
      {
        // Decode the value
        std::stringstream sword(word);
        sword >> val;
      }
      vec.push_back(val);
    }
  }

  if (nvalues != (int) vec.size())
  {
    messerr("Reading was expecting %d terms. %d found", nvalues,(int) vec.size());
    vec.clear();
    return false;
  }
  return true;
}
