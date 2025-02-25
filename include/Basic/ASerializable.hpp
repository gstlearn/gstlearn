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

#include "Basic/SerializeNeutralFile.hpp"
#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Basic/AStringable.hpp"

#include <iostream>
#include <stdarg.h>
#include <fstream>

namespace H5
{
  class Group;
};

class GSTLEARN_EXPORT ASerializable
{
public:
  ASerializable();
  ASerializable(const ASerializable& r);
  ASerializable& operator=(const ASerializable& r);
  virtual ~ASerializable();

  bool deserialize(std::istream& is, bool verbose = true);
  bool serialize(std::ostream& os,bool verbose = true) const;
  bool dumpToNF(const String& neutralFilename, bool verbose = false) const;
  bool dumpToH5(const String& H5Filename, bool verbose = false) const;

  static String buildFileName(int status, const String& filename, bool ensureDirExist = false);

  static String getHomeDirectory(const String& sub = "");
  static String getWorkingDirectory();
  static String getTestData(const String& subdir, const String& filename);
  static String getFileIdentity(const String& filename, bool verbose = false);
  static void setContainerName(bool useDefault,
                               const String& containerName = "",
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
  virtual String _getNFName() const = 0;

protected:
  virtual bool _deserialize(std::istream& is, bool verbose = false) = 0;
  virtual bool _deserializeH5(H5::Group& /*grp*/, bool /*verbose*/ = false)
  {
    // TODO virtual pure
    messerr("Not implemented yet");
    return false;
  }
  virtual bool _serialize(std::ostream& os, bool verbose = false) const = 0;
  virtual bool _serializeH5(H5::Group& /*grp*/, bool /*verbose*/ = false) const
  {
    // TODO virtual pure
    messerr("Not implemented yet");
    return false;
  }

  bool _fileOpenWrite(const String& filename,
                      std::ofstream& os,
                      bool verbose = false) const;
  bool _fileOpenRead(const String& filename,
                     std::ifstream& is,
                     bool verbose = false) const;

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
                             VectorT<T>& vec,
                             int nvalues);

  template<typename T>
  static bool _recordReadVecInPlace(std::istream& is,
                             const String& title,
                             VectorDouble::iterator& it,
                             int nvalues);

  static bool _tableRead(std::istream &is,
                         const String &string,
                         int ntab,
                         double *tab);
  static bool _tableWrite(std::ostream &os,
                          const String &string,
                          int ntab,
                          const VectorDouble &tab);

private:
  static String _myContainerName;
  static String _myPrefixName;
};

template<typename T>
bool ASerializable::_recordWrite(std::ostream& os, const String& title, const T& val)
{
  return SerializeNeutralFile::recordWrite(os, title, val);
}

template<typename T>
bool ASerializable::_recordWriteVec(std::ostream& os,
                                    const String& title,
                                    const std::vector<T>& vec)
{
  return SerializeNeutralFile::recordWriteVec(os, title, vec);
}

template<typename T>
bool ASerializable::_recordRead(std::istream& is, const String& title, T& val)
{
  return SerializeNeutralFile::recordRead(is, title, val);
}

template<typename T>
bool ASerializable::_recordReadVec(std::istream& is,
                                   const String& title,
                                   VectorT<T>& vec,
                                   int nvalues)
{
  return SerializeNeutralFile::recordReadVec(is, title, vec, nvalues);
}

template<typename T>
bool ASerializable::_recordReadVecInPlace(std::istream& is,
                                          const String& title,
                                          VectorDouble::iterator& it,
                                          int nvalues)
{
  return SerializeNeutralFile::recordReadVecInPlace<T>(is, title, it, nvalues);
}
