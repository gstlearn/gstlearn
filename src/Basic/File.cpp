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

FILE* gslFopen(const char* path, const char* mode)
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

FILE* gslFopen(const String& path, const String& mode)
{
  return gslFopen(path.c_str(), mode.c_str());
}

bool gslFileExist(const char *path, const char* mode)
{
  FILE* file = gslFopen(path, mode);
  bool exists = file != nullptr;
  if (exists) fclose(file);
  return exists;
}

bool gslFileExist(const String& path, const String& mode)
{
  FILE* file = gslFopen(path, mode);
  bool exists = file != nullptr;
  if (exists) fclose(file);
  return exists;
}

