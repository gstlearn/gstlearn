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
#include "Basic/File.hpp"
#include "Basic/String.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#ifdef __linux__ // Not operational under MacOS
#include <wordexp.h>
#endif

//#include <boost/filesystem.hpp>

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

String ASerializable::myContainerName = String();
String ASerializable::myPrefixName = String();

ASerializable::ASerializable()
{
}
/**
 * Copy constructor: don't copy temporary file info
 */
ASerializable::ASerializable(const ASerializable& /*r*/)
{
}
/**
 * Assignment operator: don't copy temporary file info
 */
ASerializable& ASerializable::operator=(const ASerializable& /*r*/)
{
  return *this;
}

ASerializable::~ASerializable()
{
}

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
  if (_fileOpenWrite(neutralFilename, os, true))
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

bool ASerializable::_fileOpenWrite(const String& filename,
                                   std::ofstream& os,
                                   bool verbose) const
{
  // Close the stream if opened
  if (os.is_open()) os.close();
  // Build the multi-platform filename
  String filepath = buildFileName(2, filename, true);
  // Open new stream
  os.open(filepath, std::ios::out | std::ios::trunc);
  if (!os.is_open())
  {
    if (verbose) messerr("Error while opening %s", filepath.c_str());
    return false;
  }
  // Write the file type (class name)
  os << _getNFName() << std::endl;
  return os.good();
}

bool ASerializable::_fileOpenRead(const String& filename,
                                  std::ifstream& is,
                                  bool verbose) const
{
  // Close the stream if opened
  if (is.is_open()) is.close();
  // Build the multi-platform filename
  String filepath = buildFileName(1, filename, true);
  // Open new stream
  is.open(filepath, std::ios::in);
  if (!is.is_open())
  {
    if (verbose) messerr("Error while opening %s", filepath.c_str());
    return false;
  }
  // Read and check the file type (class name)
  String type;
  is >> type;
  if (type != _getNFName())
  {
    if (verbose)
      messerr("The file %s has the wrong type (read: %s, expected: %s)",
              filepath.c_str(), type.c_str(), _getNFName().c_str());
    is.close();
    return false;
  }
  return is.good(); // Cannot be "end of file" already
}

bool ASerializable::_commentWrite(std::ostream& os, const String& comment)
{
  if (os.good())
  {
    if (comment.empty())
      os << std::endl;
    else
      os << "# " << comment << std::endl;
  }
  return os.good();
}

bool ASerializable::_tableWrite(std::ostream& os,
                                const String& string,
                                int ntab,
                                const VectorDouble& tab)
{
  bool ret = true;
  VectorDouble loctab(ntab);
  for (int i = 0; i < ntab; i++) loctab[i] = tab[i];
  ret = ret && _recordWriteVec<double>(os, string, loctab);
  return ret;
}

bool ASerializable::_tableRead(std::istream &is,
                               const String &string,
                               int ntab,
                               double *tab)
{
  bool ret = true;
  VectorDouble loctab(ntab);
  ret = ret && _recordReadVec<double>(is, string, loctab, ntab);
  if (!ret) return 1;
  for (int i = 0; i < ntab; i++) tab[i] = loctab[i];
  return ret;
}

bool ASerializable::_onlyBlanks(char *string)
{
  int number = static_cast<int>(strlen(string));
  for (int i = 0; i < number; i++)
  {
    if (string[i] != ' ') return false;
  }
  return true;
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
// TODO: to be restored when boost (or c++14) is usable for gstlearn (is_absolute, path manipulation, etc.)
//  boost::filesystem::path final;
//  if (! myContainerName.empty())
//  {
//    boost::filesystem::path local(myContainerName);
//    final += local;
//  }
//  if (! myPrefixName.empty())
//  {
//    boost::filesystem::path local(myPrefixName);
//    final += local;
//  }
//  boost::filesystem::path file(filename);
//  final += file;
//  String fileLocal = final.string();

  String fileLocal;

  // In the case of Output File (2), 'filename' is appended after the 'containerName' and 'prefixName'
  // In the case of Input file (1), the process depends on the contents of 'filename':
  // - if 'filename' is absolute (starts with '/' or second character is ':'): do nothing
  // - otherwise, add the 'containerName' and 'prefixName' (if defined)
  if (status == 2 || (filename.size() > 2 && filename[0] != '/' && filename[1] != ':'))
  {
    if (!myContainerName.empty())
    {
      fileLocal += myContainerName;
      if (ensureDirExist)
      {
        (void) createDirectory(fileLocal);
      }
    }
    if (!myPrefixName.empty())
    {
      fileLocal += myPrefixName;
    }
  }
  fileLocal += filename;

  String filePath = fileLocal;

#ifdef __linux__ // Not operational under MacOS
  // Check the presence of tilde character
  wordexp_t p;
  wordexp(fileLocal.c_str(), &p, 0);

  filePath = p.we_wordv[p.we_offs];
  wordfree(&p);
#endif

  return filePath;
}

String ASerializable::getHomeDirectory(const String& sub)
{
  std::stringstream sstr;
#if defined(_WIN32) || defined(_WIN64)
  String home_drive = gslGetEnv("HOMEDRIVE");
  String home_path = gslGetEnv("HOMEPATH");
  sstr << home_drive << home_path;
#else
  String home_dir = gslGetEnv("HOME");
  sstr << home_dir;
#endif
  // TODO : Cross-platform way to build file path (use boost ?)
  if (!sub.empty()) sstr << "/" << sub;
  return sstr.str();
}

String ASerializable::getWorkingDirectory()
{
  String path;
#if defined(_WIN32) || defined(_WIN64)
  char buffer[LONG_SIZE];
  if (GetModuleFileName(NULL, buffer, LONG_SIZE) != 0)
  path = String(buffer);
#else
  char buffer[LONG_SIZE];
  if (getcwd(buffer, sizeof(buffer)) != NULL) path = String(buffer);
#endif
  return path;
}

/**
 * This method returns the absolute path to a Test Data file
 * This can only be used in non-regression test (NOT in any Python or R stand-alone script)
 *
 * @param subdir Sub directory (in doc/data folder) containing the required file
 * @param filename Name of the required data file
 *
 * @return
 */
String ASerializable::getTestData(const String& subdir, const String& filename)
{
  String path = getExecDirectory();
  //std::cout << "path=" << path << std::endl;
  // TODO : Cross-platform way to build file path (use boost ?)
  // TODO : Find a proper way to register global folders (data, docs etc...)
#if defined(_WIN32) || defined(_WIN64)
  path += "..";
  path += "\\";
  path += "..";
  path += "\\";
  path += "..";
  path += "\\";
  path += "doc";
  path += "\\";
  path += "data";
  path += "\\";
  // Concatenate with the Sub-Directory
  path += subdir;
  path += "\\";
  // Concatenate with the Filename
  path += filename;
#else
  path += "../../../doc/data/";
  // Concatenate with the Sub-Directory
  path += subdir;
  path += "/";
  // Concatenate with the Filename
  path += filename;
#endif
  return path;
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
    // Default is first set to PYGSTLEARN_DIR (if defined)
    String pygst(gslGetEnv("PYGSTLEARN_DIR"));
    if (pygst.empty())
    {
      // Otherwise, it is set to HOME/gstlearn_dir
      pygst = ASerializable::getHomeDirectory("gstlearn_dir/");
      if (verbose) message("Results are stored in %s\n", pygst.c_str());
    }
    else
    {
      if (verbose)
        message("Results are stored in PYGSTLEARN_DIR\n");
    }
    myContainerName = pygst;
  }
  else
  {
    myContainerName = containerName;
  }
}

/**
 * This enables un-defining the Container Name. Then files will be saved on current Directory
 */
void ASerializable::unsetContainerName()
{
  myContainerName.erase();
}

void ASerializable::setPrefixName(const String& prefixName)
{
  myPrefixName = prefixName;
}

void ASerializable::unsetPrefixName(void)
{
  myPrefixName.erase();
}

const String& ASerializable::getContainerName()
{
  return myContainerName;
}

const String& ASerializable::getPrefixName()
{
  return myPrefixName;
}

/*!
 * Cross platform way to create a directory
 * (or ensure its existence)
 */
bool ASerializable::createDirectory(const String& dir)
{
  // TODO boost::filesystem::create_directory(dir);
#if defined(_WIN32) || defined(_WIN64)
  if (CreateDirectory(dir.c_str(), NULL) ||       // Directory creation
      ERROR_ALREADY_EXISTS == GetLastError())
  {   // or Directory was existing
    return true;
  }
  return false;
#else
  struct stat sb;
  return ((stat(dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) || // Directory exists
      (mkdir(dir.c_str(), 0755) == 0));                           // or Creation
#endif
}

/*!
 * Cross platform way to get executable directory.
 * Returned directory contains trailing separator
 */
String ASerializable::getExecDirectory()
{
  // TODO boost::filesystem::path program_location
  String dir = getHomeDirectory();
#if defined(_WIN32) || defined(_WIN64)
  char buffer[MAX_PATH] = "";
  if (GetModuleFileName(NULL, buffer, MAX_PATH) != 0)
    dir = String(buffer);
#elif defined (__APPLE__)
  char buffer[PATH_MAX] = "";
  uint32_t bufsize = PATH_MAX;
  if(!_NSGetExecutablePath(buffer, &bufsize))
    dir = String(buffer);
#else // __linux__
  char buffer[LONG_SIZE] = "";
  if (readlink("/proc/self/exe", buffer, LONG_SIZE) != -1)
    dir = String(buffer);
#endif
  return getDirectory(dir);
}

/**
 * Cross-platform way to get parent directory from a path.
 * Returned directory contains trailing separator.
 */
String ASerializable::getDirectory(const String& path)
{
  // TODO boost::filesystem::parent_path
  size_t found = path.find_last_of("/\\");
  String dir = path.substr(0, found + 1);
  return dir;
}
