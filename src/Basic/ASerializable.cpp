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
#include "Basic/AStringable.hpp"
#include "Basic/SerializeHDF5.hpp"
#include "Basic/SerializeNeutralFile.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"

#include <iostream>
#include <filesystem>
#include <fstream>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h> // for CreateDirectory
#else
#include <unistd.h> // for readlink
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h> // for _NSGetExecutablePath
#endif

#include <sys/stat.h>
#include <sys/types.h>

String ASerializable::_myContainerName = String();
String ASerializable::_myPrefixName = String();

ASerializable::ASerializable()                                    = default;
ASerializable::ASerializable(const ASerializable&)                = default;
ASerializable& ASerializable::operator=(const ASerializable&)     = default;
ASerializable::ASerializable(ASerializable&&) noexcept            = default;
ASerializable& ASerializable::operator=(ASerializable&&) noexcept = default;
ASerializable::~ASerializable()                                   = default;

bool ASerializable::deserialize(std::istream& is, bool verbose)
{
  bool ret = _deserialize(is, verbose);
  if (verbose && !ret) messerr("Problem when reading the Neutral File.");
  return ret;
}

bool ASerializable::serialize(std::ostream& os, bool verbose) const
{
  return _serialize(os, verbose);
}

bool ASerializable::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  bool ret = true;
  if (SerializeNeutralFile::fileOpenWrite(*this, neutralFilename, os, true))
  {
    ret = _serialize(os, verbose);
    if (! ret)
    {
      messerr("Problem writing in the Neutral File.");
    }
    os.close();
  }
  return ret;
}

#ifdef HDF5
bool ASerializable::dumpToH5(const String& H5Filename, bool verbose) const
{
  auto file = SerializeHDF5::fileOpenWrite(H5Filename);
  bool ret  = _serializeH5(file, verbose);
  if (!ret)
  {
    messerr("Problem writing in the HDF5 file.");
  }

  return ret;
}
#endif

bool ASerializable::_fileOpenWrite(const String& filename,
                                   std::ofstream& os,
                                   bool verbose) const
{
  return SerializeNeutralFile::fileOpenWrite(*this, filename, os, verbose);
}

bool ASerializable::_fileOpenRead(const String& filename,
                                  std::ifstream& is,
                                  bool verbose) const
{
  return SerializeNeutralFile::fileOpenRead(*this, filename, is, verbose);
}

bool ASerializable::_commentWrite(std::ostream& os, const String& comment)
{
  return SerializeNeutralFile::commentWrite(os, comment);
}

bool ASerializable::_tableWrite(std::ostream& os,
                                const String& string,
                                int ntab,
                                const VectorDouble& tab)
{
  return SerializeNeutralFile::tableWrite(os, string, ntab, tab);
}

bool ASerializable::_tableRead(std::istream &is,
                               const String &string,
                               int ntab,
                               double *tab)
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
String ASerializable::buildFileName(int status, const String& filename, bool ensureDirExist)
{
  // In the case of Output File (2), 'filename' is appended after the 'containerName' and 'prefixName'
  // In the case of Input file (1), the process depends on the contents of 'filename':
  // - if 'filename' is absolute (starts with '/' or second character is ':'): do nothing
  // - otherwise, add the 'containerName' and 'prefixName' (if defined)

  std::filesystem::path fileLocal {filename};

  if (status == 1 && fileLocal.is_absolute())
  {
    return fileLocal.string();
  }

  fileLocal.clear();

  if (!_myContainerName.empty())
  {
    fileLocal = _myContainerName;
    if (ensureDirExist)
    {
      std::filesystem::create_directory(fileLocal);
    }
  }
  const auto fname = _myPrefixName + filename;

  return (fileLocal / fname).string();
}

String ASerializable::getHomeDirectory(const String& sub)
{
#if defined(_WIN32) || defined(_WIN64)
  String home_dir = gslGetEnv("USERPROFILE");
#else
  String home_dir = gslGetEnv("HOME");
#endif
  std::filesystem::path p {home_dir};
  if (!sub.empty()) p /= sub;
  return p.string();
}

/**
 * Returns the Identity of a Neutral File which allows knowing its type
 * @param filename Name of the Neutral File
 * @param verbose Verbose flag
 * @return
 */
String ASerializable::getFileIdentity(const String& filename, bool verbose)
{
  // Preliminary check
  if (filename.empty())
  {
    if (verbose) messerr("The Neutral File Name cannot be left empty");
    return String();
  }

  // Open the File
  std::ifstream file(filename);
  if (!file.is_open())
  {
    if (verbose) messerr("Could not open the Neutral File %s", filename.c_str());
    return String();
  }

  // Read the File Header
  String filetype;
  //std::getline(file, filetype);
  gslSafeGetline(file, filetype);

  // Suppress trailing blanks
  filetype = trimRight(filetype);

  // Close the file
  file.clear();

  return filetype;
}

/**
 * Set the Container Directory Name (do not forget trailing separator "/")
 * @param useDefault True if the user wants to use automated ContainerName
 *        - defined with the global variable PYGTSLEARN_DIR
 *        - or using HOME/gstlearn_dir
 * @param containerName Name or "" for current location
 * @param verbose Verbose flag
 */
void ASerializable::setContainerName(bool useDefault,
                                     const String& containerName,
                                     bool verbose)
{
  if (useDefault)
  {
    // Default is first set to GSTLEARN_OUTPUT_DIR (if defined)
    String pygst(gslGetEnv("GSTLEARN_OUTPUT_DIR"));
    if (pygst.empty())
    {
      // Otherwise, it is set to HOME/gstlearn_dir
      pygst = ASerializable::getHomeDirectory("gstlearn_dir/");
      if (verbose) message("Results are stored in %s\n", pygst.c_str());
    }
    else
    {
      if (verbose)
        message("Results are stored in GSTLEARN_OUTPUT_DIR\n");
    }
    _myContainerName = pygst;
  }
  else
  {
    _myContainerName = containerName;
  }
}

/**
 * This enables un-defining the Container Name. Then files will be saved on current Directory
 */
void ASerializable::unsetContainerName()
{
  _myContainerName.erase();
}

void ASerializable::setPrefixName(const String& prefixName)
{
  _myPrefixName = prefixName;
}

void ASerializable::unsetPrefixName(void)
{
  _myPrefixName.clear();
}

const String& ASerializable::getContainerName()
{
  return _myContainerName;
}

const String& ASerializable::getPrefixName()
{
  return _myPrefixName;
}
