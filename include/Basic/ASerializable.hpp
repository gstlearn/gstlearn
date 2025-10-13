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

#include "Basic/AStringable.hpp"
#include "Basic/SerializeNeutralFile.hpp"
#include "Enum/EFormatNF.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <cstdarg>
#include <fstream>
#include <iostream>

namespace H5
{
class Group;
};

namespace gstlrn
{

class GSTLEARN_EXPORT ASerializable
{
public:
  ASerializable();
  ASerializable(const ASerializable& r);
  ASerializable& operator=(const ASerializable& r);
  ASerializable(ASerializable&& r) noexcept;
  ASerializable& operator=(ASerializable&& r) noexcept;
  virtual ~ASerializable();

  bool dumpToNF(const String& NFFilename,
                const EFormatNF& format = EFormatNF::fromKey("DEFAULT"),
                bool verbose            = false) const;

  static String buildFileName(Id status, const String& filename, bool ensureDirExist = false);

  static String getFileIdentity(const String& filename, bool verbose = false);
  static void setPrefixName(const String& prefixName);
  static void unsetPrefixName();
  static const String& getPrefixName();
  void setDefaultFormatNF(const EFormatNF& format);

  virtual String _getNFName() const = 0;

protected:
  virtual bool _deserializeAscii(std::istream& is, bool verbose = false) = 0;
  virtual bool _deserializeH5(H5::Group& /*grp*/, bool /*verbose*/ = false)
  {
    // TODO virtual pure
    messerr("Not implemented yet");
    return false;
  }
  virtual bool _serializeAscii(std::ostream& os, bool verbose = false) const = 0;
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
  bool _fileOpenAndDeserialize(const String& filename, bool verbose = true);

  static bool _commentWrite(std::ostream& os,
                            const String& comment);
  template<typename T>
  static bool _recordWrite(std::ostream& os,
                           const String& title,
                           const T& val);
  template<typename T>
  static bool _recordWriteVec(std::ostream& os,
                              const String& title,
                              const std::vector<T>& vec);

  template<typename T>
  static bool _recordRead(std::istream& is,
                          const String& title,
                          T& val);
  template<typename T>
  static bool _recordReadVec(std::istream& is,
                             const String& title,
                             VectorT<T>& vec,
                             Id nvalues);

  template<typename T>
  static bool _recordReadVecInPlace(std::istream& is,
                                    const String& title,
                                    VectorDouble::iterator& it,
                                    Id nvalues);

  static bool _tableRead(std::istream& is,
                         const String& string,
                         Id ntab,
                         double* tab);
  static bool _tableWrite(std::ostream& os,
                          const String& string,
                          Id ntab,
                          const VectorDouble& tab);

private:
  static String _myPrefixName;
  EFormatNF _defaultFormatNF {EFormatNF::H5};
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
                                   Id nvalues)
{
  return SerializeNeutralFile::recordReadVec(is, title, vec, nvalues);
}

template<typename T>
bool ASerializable::_recordReadVecInPlace(std::istream& is,
                                          const String& title,
                                          VectorDouble::iterator& it,
                                          Id nvalues)
{
  return SerializeNeutralFile::recordReadVecInPlace<T>(is, title, it, nvalues);
}
} // namespace gstlrn