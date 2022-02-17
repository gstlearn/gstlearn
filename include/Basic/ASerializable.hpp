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

  static String buildFileName(const String& filename, bool ensureDirExist = false);

  static String getHomeDirectory(const std::string& sub = "");
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
  virtual int _deserialize2(std::istream& s, bool verbose = false) = 0;
  virtual int _serialize2(std::ostream& s,bool verbose = false) const = 0;

  static bool _fileOpenWrite2(const String& filename,
                              const String& filetype,
                              std::ofstream& os,
                              bool verbose = false);
  static bool _fileOpenRead2(const String& filename,
                             const String& filetype,
                             std::ifstream& is,
                             bool verbose = false);

  static bool _commentWrite2(std::ostream& os,
                             const String& comment);
  template <class T>
  static bool _recordWrite2(std::ostream& os,
                            const String& title,
                            const T& val);
  template <class T>
  static bool _recordWriteVec2(std::ostream& os,
                               const String& title,
                               const T& vec);

  template <class T>
  static bool _recordRead2(std::istream& is,
                           const String& title,
                           T& val);
  template <class T>
  static bool _recordReadVec2(std::istream& is,
                              const String& title,
                              T& vec);


  //virtual int _deserialize(FILE* file, bool verbose = false) = 0;
  //virtual int _serialize(FILE* file, bool verbose = false) const = 0;

  static FILE* _fileOpen(const String& filename,
                         const String& filetype,
                         const String& mode,
                         bool verbose = false);
  static int _fileClose(FILE* file, bool verbose = false);
  static int _recordRead(FILE* file,
                         const char* title,
                         const char* format, ...);
  static void _recordWrite(FILE* file, const char* format, ...);
  static void _tableWrite(FILE *file, const String& string, int ntab, const double *tab);
  static int  _tableRead(FILE* file, int ntab, double *tab);
  static int  _fileRead(FILE* file, const String& format, va_list ap);
  static void _fileWrite(FILE* file, const String& format, va_list ap);
  static bool _onlyBlanks(char *string);

  static int  _tableRead2(std::istream& is, int ntab, double *tab);
  static bool _tableWrite2(std::ostream& os,
                           const String& string,
                           int ntab,
                           const VectorDouble& tab);

private:
  static String myContainerName;
  static String myPrefixName;
};

template <class T>
bool ASerializable::_recordWrite2(std::ostream& os,
                                  const String& title,
                                  const T& val)
{
  if (os.good())
  {
    if (isNA(val))
    {
      if (title.empty())
         os << "NA" << " ";
       else
         os << "NA" << " # " << title << std::endl;
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

template <class T>
bool ASerializable::_recordWriteVec2(std::ostream& os,
                                     const String& title,
                                     const T& vec)
{
  if (os.good())
  {
    if (!title.empty())
      os << "# " << title << std::endl;
    for (auto val : vec)
    {
//      if (isNA(&val))
//        os << "NA" << " ";
//      else
        os << val << " ";
    }
    os << std::endl;
  }
  return os.good();
}

template <class T>
bool ASerializable::_recordRead2(std::istream& is,
                                 const String& title,
                                 T& val)
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
      word = trimLeft(word);
      if (!word.empty())
      {
        if (word[0] != '#')
          break;   // We found something
        else
          std::getline(is, word); // Eat all the comment
      }
    }

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
  return is.good();
}

template <class T>
bool ASerializable::_recordReadVec2(std::istream& is,
                                    const String& title,
                                    T& vec)
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
      line = trimLeft(line);
      if (!line.empty() && line[0] != '#')
        break;
    }
    // Decode the line
    std::stringstream sstr(line);
    while (sstr.good())
    {
      typename T::value_type val;
      sstr >> val;
      if (!sstr.good() && !sstr.eof())
      {
        messerr("Error while reading %s", title.c_str());
        vec.clear();
        return false;
      }
      if (sstr.good()) vec.push_back(val);
    }
  }
  return is.good();
}
