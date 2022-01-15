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

#include <stdarg.h>
#include <stdio.h>

class GSTLEARN_EXPORT ASerializable
{
public:
  ASerializable();
  ASerializable(const ASerializable& r);
  ASerializable& operator=(const ASerializable& r);
  virtual ~ASerializable();

  static String buildFileName(const String& filename, bool ensureDirExist = false);

  static String getHomeDirectory(const std::string& sub = "");
  static String getFileIdentify(const String& filename);
  static void setContainerName(bool useDefault,
                               const String& containerName = String(),
                               bool verbose = false);
  static void setPrefixName(const String& prefixName);
  static const String& getContainerName();
  static const String& getPrefixName();

  // TODO : Directory manipulation class
  static bool createDirectory(const String& dir);
  static String getExecDirectory();
  static String getDirectory(const String& path);

protected:
  virtual int _deserialize(FILE* file, bool verbose = false) = 0;
  virtual int _serialize(FILE* file, bool verbose = false) const = 0;

  static FILE* _fileOpen(const String& filename,
                         const String& filetype,
                         const String& mode,
                         bool verbose = false);
  static int _fileClose(FILE* file, bool verbose = false);
  static int _recordRead(FILE* file,
                         const char* title,
                         const char* format, ...);
  static void _recordWrite(FILE* file, const char* format, ...);
  static int _fileRead(FILE* file, const String& format, va_list ap);
  static void _fileWrite(FILE* file, const String& format, va_list ap);
  static bool _onlyBlanks(char *string);

private:
  static String myContainerName;
  static String myPrefixName;
};
