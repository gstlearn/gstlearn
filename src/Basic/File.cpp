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

char* gslStrcpy(char* dst, int dst_size, const char* src)
{
  char *res;
#ifdef __unix
  res = strcpy(dst, src);
#else
  int err;
  err = strcpy_s(dst, dst_size, src);
  res = dst;
#endif
  return res;
}
void gslStrcpy(String& dst, const String& src)
{
  dst = src;
}

char* gslStrcat(char* dst, int dst_size, const char* src)
{
  char *res;
#ifdef __unix
  res = strcat(dst, src);
#else
  int err;
  err = strcat_s(dst, dst_size, src);
  res = dst;
#endif
  return res;
}

void gslStrcat(String& dst, const String& src)
{
  dst = dst + src;
}

int gslSPrintf(char* dst, int dst_size, const char* format, ...)
{
  int ret;
  va_list ap;
  va_start(ap,format);

#if __STDC_WANT_LIB_EXT1__ == 1
  ret = vsprintf_s( dst, dst_size, format, ap);
#else
  ret = vsprintf( dst, format, ap);
#endif
  va_end(ap);
  return ret;
}

int gslSPrintf(String& dst, const String& format, ...)
{
  int size=100;
  va_list ap;

  while (1) {
      dst.resize(size);
      va_start(ap, format);
      int n = vsnprintf(&dst[0], size, format.c_str(), ap);
      va_end(ap);

      if (n > -1 && n < size) {
          dst.resize(n); // Make sure there are no trailing zero char
          return n;
      }
      if (n > -1)
          size = n + 1;
      else
          size *= 2;
  }
}
