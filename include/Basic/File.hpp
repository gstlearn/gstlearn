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

/**
 * Redirection facility
 * This function is used to redirect the output of the script to an auxiliary file
 * This facility is used within the non-regression tests calling:
 *       std::stringstream sfn;
 *       sfn << gslBaseName(__FILE__) << ".out";
 *       StdoutRedirect sr(sfn.str(), argc, argv);
 * The redirection can be cancelled (with no argument modification)
 * if the number of arguments 'argc' is larger than a given threshold
 * By default, this threshold is set to 1 (name of executable itself)
 * This can even be modified by adding 'number' different from 1 if some arguments
 * are compulsory.
 */
class GSTLEARN_EXPORT StdoutRedirect {
public:
  StdoutRedirect(const String &file = "",
                 int argc = 0,
                 char *argv[] = nullptr,
                 int number = 1);
  ~StdoutRedirect();
  StdoutRedirect(const StdoutRedirect&) = delete;
  StdoutRedirect& operator=(const StdoutRedirect&) = delete;

  void start(const String& file);
  void stop();

private:
  bool _flagActive;
#if defined(_WIN32) || defined(_WIN64)
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
GSTLEARN_EXPORT void skipBOM(std::ifstream &ins);
GSTLEARN_EXPORT FILE* gslFopen(const char *path, const char* mode);
GSTLEARN_EXPORT FILE* gslFopen(const String& path, const String& mode);
GSTLEARN_EXPORT bool gslFileExist(const char *path, const char* mode);
GSTLEARN_EXPORT bool gslFileExist(const String& path, const String& mode);

GSTLEARN_EXPORT String gslBaseName(const String& path, bool keepExtension = false);

GSTLEARN_EXPORT String gslGetEnv(const String& name);

/**
 * Get line from an input text stream whatever the end of line convention.
 * Thanks to: https://stackoverflow.com/a/6089413
 */
GSTLEARN_EXPORT std::istream& gslSafeGetline(std::istream& is, String& t);

