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

#include <fstream>

/*
// TODO: If you restore that (for HANDLE type below) you obtain an
//       undefined reference of EDbg::IPrintDialogServices in io.cpp !!!!
//       Does 'EDbg' short name is conflicting under Windows?
#if defined(_WIN32) || defined(_WIN64)
#  include <windows.h>
#endif
*/

// Standard output stream redirection

class GSTLEARN_EXPORT StdoutRedirect {
public:
  StdoutRedirect(const String& file = "");
  ~StdoutRedirect();
  StdoutRedirect(const StdoutRedirect&) = delete;
  StdoutRedirect& operator=(const StdoutRedirect&) = delete;

  void start(const String& file);
  void stop();

private:
#ifdef _WIN32
  // HANDLE _old_stdout;
  void* _old_stdout;
#else
  std::streambuf* _coutbuf;
  std::ofstream _out;
#endif
};

/**
 * This file contains all the OLD-STYLE declarations causing warnings on Windows
 * They should gradually be replaced by modern statements
 * However, they are kept there to keep track on these statements
 */

// TODO : File manipulation class

// Skips the Byte Order Mark (BOM) that defines UTF-8 in some text files.
//https://stackoverflow.com/a/17219495
GSTLEARN_EXPORT void skipBOM(std::ifstream &in);
GSTLEARN_EXPORT FILE* gslFopen(const char *path, const char* mode);
GSTLEARN_EXPORT FILE* gslFopen(const String& path, const String& mode);
GSTLEARN_EXPORT bool gslFileExist(const char *path, const char* mode);
GSTLEARN_EXPORT bool gslFileExist(const String& path, const String& mode);

GSTLEARN_EXPORT String gslBaseName(const String& path, bool keepExtension = false);

GSTLEARN_EXPORT String gslGetEnv(const String& name);

