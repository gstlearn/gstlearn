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
#include "Basic/File.hpp"
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdarg.h>
#include <stdlib.h>

// Skips the Byte Order Mark (BOM) that defines UTF-8 in some text files.
//https://stackoverflow.com/a/17219495
GSTLEARN_EXPORT void skipBOM(std::ifstream &in)
{
  char test[3] = {0};
  in.read(test, 3);
  if ((unsigned char)test[0] == 0xEF &&
      (unsigned char)test[1] == 0xBB &&
      (unsigned char)test[2] == 0xBF)
  {
      return;
  }
  in.seekg(0);
}

GSTLEARN_EXPORT FILE* gslFopen(const char* path, const char* mode)
{
  FILE *file;

#ifdef __unix
  file = fopen(path, mode);
#else
  errno_t err;
  err = fopen_s( &file, path, mode );
  if ( err != 0 ) return nullptr;
#endif
  return file;
}

GSTLEARN_EXPORT FILE* gslFopen(const String& path, const String& mode)
{
  return gslFopen(path.c_str(), mode.c_str());
}

GSTLEARN_EXPORT bool gslFileExist(const char *path, const char* mode)
{
  FILE* file = gslFopen(path, mode);
  bool exists = file != nullptr;
  if (exists) fclose(file);
  return exists;
}

GSTLEARN_EXPORT bool gslFileExist(const String& path, const String& mode)
{
  FILE* file = gslFopen(path, mode);
  bool exists = file != nullptr;
  if (exists) fclose(file);
  return exists;
}

GSTLEARN_EXPORT char* gslGetEnv(const char* name)
{
#if defined(_WIN32) || defined(_WIN64)
  return getenv(name);
#else
  return std::getenv(name);
#endif
}

