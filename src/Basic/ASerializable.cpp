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
#include "Basic/ASerializable.hpp"
#include "Basic/File.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/SerializeNeutralFile.hpp"
#include "Basic/String.hpp"
#include "Enum/EFormatNF.hpp"

#include <clocale>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace gstlrn
{
  String ASerializable::_myPrefixName = String();
  thread_local bool locale_set = false;

  ASerializable::ASerializable()
  {
    if (!locale_set)
    {
      // support utf8 characters in Windows files/folders
      setlocale(LC_ALL, "en_US.utf8");
      locale_set = true;
    }
  }

  void ASerializable::setDefaultFormatNF(const EFormatNF& format)
  {
    _defaultFormatNF = format;
  }

  /**
   * @brief Dump the contents of an object into an Output File
   * using a given Output NF Format
   *
   * @param NFFilename Name of the Output File
   * @param format Choice of the format (see remarks)
   * @return true or false
   *
   * @remarks In the argument 'format', the user can select the format for encoding
   * the contents of the output file.
   * If the value DEFAULT is used, the package uses the Format currently defined
   * as the defaulted one. This default value can be updated using the method
   * ASerializable::DefaultFormatNF()
   */
  bool ASerializable::dumpToNF(
    const String& NFFilename,
    const EFormatNF& format) const
  {
    bool ret = true;

    EFormatNF formatLocal = format;
    if (format == EFormatNF::DEFAULT)
    {
      formatLocal = _defaultFormatNF;
    }

    // Check if H5 format is available: otherwise force ASCII
#ifndef HDF5
    formatLocal = EFormatNF::ASCII;
#endif

    if (formatLocal == EFormatNF::ASCII)
    {
      std::ofstream os;
      if (SerializeNeutralFile::fileOpenWrite(*this, NFFilename, os, true))
      {
        ret = _serializeAscii(os);
        if (!ret) messerr("Problem writing in the Neutral File.");
        os.close();
      }
      return ret;
    }

#ifdef HDF5
    if (formatLocal == EFormatNF::H5)
    {
      auto file = SerializeHDF5::fileOpenWrite(*this, NFFilename);
      ret = serializeH5(file);
      if (!ret) messerr("Problem writing in the HDF5 file.");
      return ret;
    }
#endif

    messerr(
      "Aserializable::dumpToNF: No Format is defined (enum value: %d)",
      formatLocal.getValue());
    return false;
  }

  bool
    ASerializable::_fileOpenAndDeserialize(const String& filename, bool verbose)
  {
    // Check that the file exists
    String filepath = ASerializable::buildFileName(1, filename, true);
    std::ifstream file(filepath);
    if (!file.good())
    {
      if (verbose) messerr("The file %s does not exist", filepath.c_str());
      return false;
    }

    // Try to open it according to HDF5 format
#ifdef HDF5
    if (H5::H5File::isHdf5(filepath))
    {
      auto fileStr = SerializeHDF5::fileOpenRead(filename);

      return deserializeH5(fileStr);
    }
#endif

    // Try to open it according to ASCII format
    std::ifstream is;
    if (SerializeNeutralFile::fileOpenRead(*this, filename, is, verbose))
    {
      return _deserializeAscii(is);
    }

    if (verbose) messerr("Opening the file %s failed", filename.c_str());
    return false;
  }

  bool ASerializable::_fileOpenWrite(
    const String& filename,
    std::ofstream& os,
    bool verbose) const
  {
    return SerializeNeutralFile::fileOpenWrite(*this, filename, os, verbose);
  }

  bool ASerializable::_fileOpenRead(
    const String& filename,
    std::ifstream& is,
    bool verbose) const
  {
    return SerializeNeutralFile::fileOpenRead(*this, filename, is, verbose);
  }

  bool ASerializable::_commentWrite(std::ostream& os, const String& comment)
  {
    return SerializeNeutralFile::commentWrite(os, comment);
  }

  bool ASerializable::_tableWrite(
    std::ostream& os,
    const String& string,
    Id ntab,
    const VectorDouble& tab)
  {
    return SerializeNeutralFile::tableWrite(os, string, ntab, tab);
  }

  bool ASerializable::_tableRead(
    std::istream& is,
    const String& string,
    Id ntab,
    double* tab)
  {
    return SerializeNeutralFile::tableRead(is, string, ntab, tab);
  }

  /**
   * Build a standard filename for Read or Write operation
   * @param status 1 for Read and 2 for Write
   * @param filename Name of the filename (see remark)
   * @param ensureDirExist When TRUE, the Directory is created if not already existing
   * @return
   */
  String ASerializable::buildFileName(
    Id status,
    const String& filename,
    bool ensureDirExist)
  {
    // In the case of Output File (2), 'filename' is appended after the 'containerName' and 'prefixName'
    // In the case of Input file (1), the process depends on the contents of 'filename':
    // - if 'filename' is absolute (starts with '/' or second character is ':'): do nothing
    // - otherwise, add the 'containerName' and 'prefixName' (if defined)

    std::filesystem::path fileLocal{filename};

    if (status == 1 && fileLocal.is_absolute())
    {
      return fileLocal.string();
    }
    fileLocal.clear();

    // container name: first search for the GSTLEARN_OUTPUT_DIR
    // environment variable. If not defined, use the current directory
    // instead.
    const auto output_dir = gslGetEnv("GSTLEARN_OUTPUT_DIR");

    if (!output_dir.empty())
    {
      fileLocal = output_dir;
    }
    else
    {
      fileLocal = std::filesystem::current_path();
    }

    if (ensureDirExist)
    {
      std::filesystem::create_directory(fileLocal);
    }

    const auto fname = _myPrefixName + filename;

    return (fileLocal / fname).string();
  }

  /**
   * Returns the Identity of a Neutral File which allows knowing its type
   * @param filename Name of the Neutral File
   * @param verbose Verbose flag
   * @return
   */
  String ASerializable::getFileIdentity(const String& filename, bool verbose)
  {
    // Preliminary check (no message if string is empty ... even in verbose)
    if (filename.empty()) return String();
    if (verbose) message("Input File Name = %s\n", filename.c_str());

    // Build the multi-platform filename
    const auto filepath = ASerializable::buildFileName(1, filename, true);
    if (verbose) message("Input File Path = %s\n", filepath.c_str());

    // Open the file according to various formats
    Id ret_type = -1;
    String classType;

#ifdef HDF5
    if (ret_type < 0)
    {
      // Check if the file is written at format HDF5
      if (H5::H5File::isHdf5(filepath))
      {

        // Attempt to open according to a H5 format
        H5::H5File file{filepath, H5F_ACC_RDONLY};

        if (!file.nameExists("gstlearn metadata"))
        {
          messerr(
            "File %s doesn't contain Gstlearn metadata…", filepath.c_str());
          return String();
        }
        ret_type = 1;

        // Read the class type
        classType = SerializeHDF5::getFileClass(filepath, verbose);
      }
    }
#endif

    if (ret_type < 0)
    {

      // Attempt to open according to the ASCII format
      std::ifstream file(filepath);
      if (!file.is_open())
      {
        if (verbose)
          messerr("Could not open the Neutral File %s", filepath.c_str());
        return String();
      }
      ret_type = 0;

      // Read the Class Type
      gslSafeGetline(file, classType);
      classType = trimRight(classType);

      // Close the file
      file.clear();
    }

    if (verbose) message("Decoded Type = %s\n", classType.c_str());

    return classType;
  }

  void ASerializable::setPrefixName(const String& prefixName)
  {
    _myPrefixName = prefixName;
  }

  void ASerializable::unsetPrefixName(void)
  {
    _myPrefixName.clear();
  }

  const String& ASerializable::getPrefixName()
  {
    return _myPrefixName;
  }

#ifdef HDF5
  /**
   * @brief Returns the Full Path of a given HDF5 Group
   *
   * @param group   HDF5 Group
   * @return String
   */
  String ASerializable::getGroupFullPath(const H5::Group& group)
  {
    ssize_t size = H5Iget_name(group.getId(), nullptr, 0);
    std::string path(size, '\0');
    H5Iget_name(group.getId(), path.data(), size + 1);
    return path;
  }

  /**
   * @brief Returns the list of Parent Groups of a given HDF5 Group
   *
   * @param group   HDF5 Group
   * @return VectorString
   */
  VectorString ASerializable::getGroupParents(const H5::Group& group)
  {
    String path = getGroupFullPath(group);
    VectorString parents;

    while (path != "/" && !path.empty())
    {
      parents.push_back(path);
      auto pos = path.find_last_of('/');
      if (pos == 0)
        path = "/";
      else
        path = path.substr(0, pos);
    }

    parents.push_back("/"); // racine
    return parents;
  }
#endif

} // namespace gstlrn
