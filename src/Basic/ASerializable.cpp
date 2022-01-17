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
#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <regex>
//#include <boost/filesystem.hpp>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h> // for CreateDirectory
#else
#include <unistd.h> // for readlink
#endif

#include <sys/stat.h>
#include <sys/types.h>

static char LINE[LONG_SIZE], LINE_MEM[LONG_SIZE], *LCUR;
static char *cur = NULL;

String ASerializable::myContainerName = String();
String ASerializable::myPrefixName    = String();

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

/****************************************************************************/
/*!
 **   Open an ASCII file
 **
 ** \return  FILE returned pointer
 **
 ** \param[in]  filename Local file name
 ** \param[in]  filetype Type of the file (optional [NULL] when 'w')
 ** \param[in]  mode     "r" or "w"
 ** \param[in]  verbose  Verbose flag
 **
 *****************************************************************************/
FILE* ASerializable::_fileOpen(const String& filename,
                               const String& filetype,
                               const String& mode,
                               bool verbose)
{
  FILE* file = nullptr;

  if (filename.empty() || filetype.empty())
  {
    messerr("Arguments 'filename' or 'filetype' are not defined");
    return nullptr;
  }
  else
  {
    // Build the multi-platform filename and open it
    String fileComplete = buildFileName(filename, true);

    // Optional printout
    if (verbose)
    {
      message("Attempt to Open the File (%s) in mode (%s)\n",
              fileComplete.c_str(), mode.c_str());
    }

    file = gslFopen(fileComplete, mode);
    if (file == nullptr)
    {
      if (verbose)
        messerr("Error when opening the Neutral File %s", fileComplete.c_str());
      return nullptr;
    }
    else
    {
      if (verbose)
        message("File successfully opened\n");
    }

    // Preliminary action
    if (mode == "r")
    {
      char idtype[LONG_SIZE];
      if (_recordRead(file, "File Type", "%s", idtype))
      {
        _fileClose(file, false);
        return nullptr;
      }
      if (strcmp(idtype,filetype.c_str()))
      {
        messerr(
            "Error: in the File (%s), its Type (%s) does not match the requested one (%s)",
            filename.c_str(), idtype, filetype.c_str());
        _fileClose(file, false);
        return nullptr;
      }
    }
    else
    {
      // Write the File ID
      if (!filetype.empty())
      {
        _recordWrite(file, "%s", filetype.c_str());
        _recordWrite(file, "\n");
      }
    }
  }
  return file;
}

int ASerializable::_fileClose(FILE *file, bool verbose)
{
  if (file != nullptr)
    fclose(file);

  // Optional printout

  if (verbose)
    message("File is successfully closed\n");

  return 0;
}

/**
 *
 * @param file  FILE structure
 * @param title Title to be printed
 * @param format String to be completed
 * @return
 *
 * @remarks: Format is not a reference here:
 * https://stackoverflow.com/questions/222195/are-there-gotchas-using-varargs-with-reference-parameters
 * => roll back to const char* to prevent warning C4840 (MSVC)
 */
int ASerializable::_recordRead(FILE* file, const char* title, const char* format, ...)
{
  va_list ap;
  int error;

  error = 0;
  va_start(ap, format);
  error = _fileRead(file, format, ap);

  if (error > 0)
    messerr("Error when reading '%s'", title);

  va_end(ap);
  return (error);
}

/**
 * Record a String in the Serialized file
 * @param file FILE structure
 * @param format String to be completed
 * @param ... Variable list of arguments
 * @return
 *
 * @remark Format is not a reference here:
 * https://stackoverflow.com/questions/222195/are-there-gotchas-using-varargs-with-reference-parameters
 * => roll back to const char* to prevent warning C4840 (MSVC)
 *
 * TODO : Impose that va_list arguments are stringable ? (For example, we can serialize ECov objects !)
 */
void ASerializable::_recordWrite(FILE* file, const char* format, ...)
{
  va_list ap;
  va_start(ap, format);
  _fileWrite(file, format, ap);
  va_end(ap);
}

/****************************************************************************/
/*!
 **  Read the next token from the file
 **
 ** \return  -1 if the end-of-file has been found
 ** \return   1 for a decoding error
 ** \return   0 otherwise
 **
 ** \param[in]  file       FILE structure
 ** \param[in]  format     format
 ** \param[in]  ap         Value to be read
 **
 ** TODO : template function
 **
 *****************************************************************************/
int ASerializable::_fileRead(FILE* file, const String& format, va_list ap)
{
  char DEL_COM = '#';
  char DEL_SEP = ' ';
  char DEL_BLK = ' ';
  char DEL_TAB = '\t';

  const char *fmt;
  int    *ret_i;
  float  *ret_f;
  double *ret_d;
  char   *ret_s;

  /* Loop on the elements to read (from the format) */

  int ideb = 0;
  while (ideb < static_cast<int>(format.size()))
  {
    /* Eliminate the blanks */

    if (format[ideb] == DEL_BLK)
    {
      ideb++;
      continue;
    }

    label_start: fmt = &format[ideb];
    if (LCUR == NULL)
    {

      /* Read the next line */

      if (fgets(LINE, LONG_SIZE, file) == NULL) return (-1);
      size_t wsize = strcspn(LINE, "\r\n");
      if (wsize > 0)
        LINE[wsize] = 0;
      else
        LINE[strlen(LINE)-1] = '\0';
      (void) gslStrcpy(LINE_MEM, LINE);

      /* Eliminate the comments and replace <TAB> by blank*/

      int flag_com = 0;
      for (unsigned int i = 0; i < strlen(LINE); i++)
      {
        if (LINE[i] == DEL_TAB) LINE[i] = DEL_BLK;
        if (LINE[i] == DEL_COM)
        {
          flag_com = 1 - flag_com;
          LINE[i] = '\0';
        }
        else
        {
          if (flag_com) LINE[i] = '\0';
        }
      }
      cur = LINE;
    }

    /* Decode the line looking for the next token */

    LCUR = gslStrtok(cur, &DEL_SEP);
    cur = NULL;
    if (LCUR == NULL) goto label_start;

    /* Reading */

    if (!strcmp(fmt, "%s"))
    {
      ret_s = va_arg(ap, char *);
      if (!_onlyBlanks(LCUR))
      {
        if (gslSScanf(LCUR, "%s", ret_s) <= 0) return (1);
      }
      ideb += 2;
    }
    else if (!strcmp(fmt, "%d"))
    {
      ret_i = va_arg(ap, int *);
      if (gslSScanf(LCUR, "%d", ret_i) <= 0) return (1);
      ideb += 2;
      if (*ret_i == (int) ASCII_TEST) *ret_i = ITEST;
    }
    else if (!strcmp(fmt, "%f"))
    {
      ret_f = va_arg(ap, float *);
      if (gslSScanf(LCUR, "%f", ret_f) <= 0) return (1);
      ideb += 2;
      if (*ret_f == ASCII_TEST) *ret_f = TEST;
    }
    else if (!strcmp(fmt, "%lf"))
    {
      ret_d = va_arg(ap, double *);
      if (gslSScanf(LCUR, "%lf", ret_d) <= 0) return (1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
    }
    else if (!strcmp(fmt, "%lg"))
    {
      ret_d = va_arg(ap, double *);
      if (gslSScanf(LCUR, "%lg", ret_d) <= 0) return (1);
      ideb += 3;
      if (*ret_d == ASCII_TEST) *ret_d = TEST;
    }
    else
    {
      messerr("Wrong format %s", fmt);
      va_end(ap);
      return (2);
    }
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Write the next token from the file
 **
 ** \param[in]  file       FILE structure
 ** \param[in]  format     Encoding format
 ** \param[in]  ap         Value to be written
 **
 ** TODO : template function
 **
 *****************************************************************************/
void ASerializable::_fileWrite(FILE* file, const String& format, va_list ap)
{
  double ret_d;
  char *ret_s;
  bool no_blank = false;

  /* Writing */

  if (format == "%s")
  {
    ret_s = va_arg(ap, char *);
    fprintf(file, "%s", ret_s);
  }
  else if (format == "%d")
  {
    int ret_i = va_arg(ap, int);
    if (FFFF(ret_i))
      fprintf(file, "%5.1lf", ASCII_TEST);
    else
      fprintf(file, "%d", ret_i);
  }
  else if (format == "%f")
  {
    ret_d = va_arg(ap, double);
    if (FFFF(ret_d))
      fprintf(file, "%5.1lf", ASCII_TEST);
    else
      fprintf(file, "%f", ret_d);
  }
  else if (format == "%lf")
  {
    ret_d = va_arg(ap, double);
    if (FFFF(ret_d))
      fprintf(file, "%5.1lf", ASCII_TEST);
    else
      fprintf(file, "%lf", ret_d);
  }
  else if (format == "%lg")
  {
    ret_d = va_arg(ap, double);
    if (FFFF(ret_d))
      fprintf(file, "%5.1lf", ASCII_TEST);
    else
      fprintf(file, "%lg", ret_d);
  }
  else if (format == "\n")
  {
    fprintf(file, "\n");
    no_blank = true;
  }
  else if (format == "#")
  {
    ret_s = va_arg(ap, char *);
    fprintf(file, "# %s\n", ret_s);
    no_blank = true;
  }
  else
  {
    messerr("Wrong format %s", format.c_str());
    return;
  }
  if (!no_blank) fprintf(file, " ");
  return;
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

String ASerializable::buildFileName(const String& filename, bool ensureDirExist)
{
// TODO: to be restored when boost is usable for pygstlearn
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
  if (!myContainerName.empty())
  {
    fileLocal += myContainerName;
    if (ensureDirExist)
    {
      (void)createDirectory(fileLocal);
    }
  }
  if (!myPrefixName.empty())
  {
    fileLocal += myPrefixName;
  }
  fileLocal += filename;

  return fileLocal;
}

String ASerializable::getHomeDirectory(const String& sub)
{
  // https://stackoverflow.com/a/2552458
#if defined(_WIN32) || defined(_WIN64)
  char* HomeDirectory = gslGetEnv("HOMEDIR");
  const char* Homepath = gslGetEnv("HOMEPATH");
  size_t size = strlen(HomeDirectory) + strlen(Homepath) + 1;
  HomeDirectory = static_cast<char *>(malloc(size));
  gslStrcat(HomeDirectory, Homepath);
#else
  const char* HomeDirectory = gslGetEnv("HOME");
#endif
  std::stringstream sstr;
  // TODO : Cross-platform way to build file path (use boost ?)
  sstr << String(HomeDirectory);
  if (!sub.empty())
    sstr << "/" << sub;
  return sstr.str();
}

/**
 * This method returns the absolute path to a Test Data
 * @return
 */
String ASerializable::getTestData(const String& subdir, const String& filename)
{
  String dirname = getExecDirectory();
#if defined(_WIN32) || defined(_WIN64)
  dirname += "\\";
  dirname += "..";
  dirname += "\\";
  dirname += "..";
  dirname += "\\";
  dirname += "doc";
  dirname += "\\";
  dirname += "data";
  dirname += "\\";
#else
  dirname += "../../doc/data/";
#endif
  std::stringstream sstr;
  // TODO : Cross-platform way to build file path (use boost ?)

  // Concatenate with the Sub-Directory (if defined)

  sstr << String(dirname);
  sstr << subdir << "/";

  // Concatenate with the Filename
  sstr << filename;

  return sstr.str();
}

String ASerializable::getFileIdentify(const String& filename)
{
  // Preliminary check
  if (filename.empty())
  {
    messerr("The Neutral File Name cannot be left empty");
    return String();
  }

  // Open the File
  std::ifstream file(filename);
  if (!file.is_open())
  {
    messerr("Could not open the Neutral File %s",filename.c_str());
    return String();
  }

  // Read the File Header
  String filetype;
  std::getline(file, filetype);

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
    char* pydir(gslGetEnv("PYGSTLEARN_DIR"));
    String pygst;
    if (pydir == nullptr)
    {
      // Otherwise, it is set to HOME/gstlearn_dir
      pygst = ASerializable::getHomeDirectory("gstlearn_dir/");
      if (verbose)
        std::cout << "Results are stored in" << pygst << std::endl;
    }
    else
    {
      pygst = pydir;
      if (verbose)
        std::cout << "Results are stored in PYGSTLEARN_DIR" << std::endl;
    }
    myContainerName = pygst;
  }
  else
  {
    myContainerName = containerName;
  }
}

void ASerializable::setPrefixName(const String& prefixName)
{
  myPrefixName = prefixName;
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
      ERROR_ALREADY_EXISTS == GetLastError()) {   // or Directory was existing
    return true;
  }
#else
  struct stat sb;
  if ((stat(dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) || // Directory exists
      (mkdir(dir.c_str(), 0755) == 0))                        // or Creation
    return true;
#endif
  return false;
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
  char buffer[MAX_PATH];
  if (GetModuleFileName(NULL, buffer, MAX_PATH) != 0)
    dir = String(buffer);
#else
  char buffer[LONG_SIZE];
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
  String dir = path.substr(0,found+1);
  return dir;
}
