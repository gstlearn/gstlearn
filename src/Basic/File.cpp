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
#include "Basic/File.hpp"
#include "Basic/AStringable.hpp"
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <iostream>

#include <fstream>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#endif

StdoutRedirect::StdoutRedirect(const String& file, int argc, char* argv[]) :
  _flagActive(true),
#if defined(_WIN32) || defined(_WIN64)
  _old_stdout(0)
#else
  _coutbuf(nullptr),
  _out()
#endif
{
  DECLARE_UNUSED(argv);
  _flagActive = (argc <= 1);
  if (!file.empty() && _flagActive)
    start(file);
}

StdoutRedirect::~StdoutRedirect()
{
  if (_flagActive)
    stop();
}

/**
 * Save current stdout handle and redirect std::cout to a file
 *
 * @param[in] file File path to be written
 */
void StdoutRedirect::start(const String& file)
{
#if defined(_WIN32) || defined(_WIN64)
  // https://stackoverflow.com/questions/54094127/redirecting-stdout-in-win32-does-not-redirect-stdout/54096218
  _old_stdout = GetStdHandle(STD_OUTPUT_HANDLE);
  HANDLE new_stdout = CreateFileA(file.c_str(), GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
  SetStdHandle(STD_OUTPUT_HANDLE, new_stdout);
  int fd = _open_osfhandle((intptr_t)new_stdout, _O_WRONLY|_O_TEXT);
  _dup2(fd, _fileno(stdout));
  _close(fd);
#else
  _coutbuf = std::cout.rdbuf();
  _out.open(file, std::fstream::out | std::fstream::trunc);
  std::cout.rdbuf(_out.rdbuf());
#endif
}

/**
 *  Restore original stdout
 */
void StdoutRedirect::stop()
{
#if defined(_WIN32) || defined(_WIN64)
  // https://stackoverflow.com/questions/32185512/output-to-console-from-a-win32-gui-application-on-windows-10
  SetStdHandle(STD_OUTPUT_HANDLE, _old_stdout);
  int fd = _open_osfhandle((intptr_t)_old_stdout, _O_WRONLY|_O_TEXT);
  if (fd >= 0) // fd could be negative for an unknown reason (https://github.com/gstlearn/gstlearn/issues/111)
  {
    FILE* fp = _fdopen(fd, "w");
    freopen_s( &fp, "CONOUT$", "w", stdout);
  }
#else
  std::cout.rdbuf(_coutbuf);
  _out.close();
#endif
}

// Skips the Byte Order Mark (BOM) that defines UTF-8 in some text files.
//https://stackoverflow.com/a/17219495
void skipBOM(std::ifstream &in)
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

FILE* gslFopen(const char* path, const char* mode)
{
  FILE *file;

#ifdef __unix
  file = fopen(path, mode);
#elif defined(__APPLE__)
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

String gslBaseName(const String& path, bool keepExtension)
{
  // TODO c++17 : use #include <filesystem> to retrieve base name
  // https://stackoverflow.com/questions/8520560/get-a-file-name-from-a-path
  String base_filename = path.substr(path.find_last_of("/\\") + 1);
  if (!keepExtension)
  {
    String::size_type p(base_filename.find_last_of('.'));
    base_filename = base_filename.substr(0, p);
  }
  return base_filename;
}


String gslGetEnv(const String& name)
{
  String text;
#if defined(_WIN32) || defined(_WIN64)
  const DWORD buffSize = 65535;
  static char buffer[buffSize];
  if (GetEnvironmentVariable(name.c_str(), buffer, buffSize))
    text = String(buffer);
#elif defined(__linux__) || defined(__APPLE__)
  char* value = std::getenv(name.c_str());
  if (value != NULL)
    text = String(value);
#endif
  return text;
}

std::istream& gslSafeGetline(std::istream& is, String& t)
{
  t.clear();

  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.

  std::istream::sentry se(is, true);
  std::streambuf *sb = is.rdbuf();

  for (;;)
  {
    int c = sb->sbumpc();
    switch (c)
    {
      case '\n':
        return is;
      case '\r':
        if (sb->sgetc() == '\n') sb->sbumpc();
        return is;
      case std::streambuf::traits_type::eof():
        // Also handle the case when the last line has no line ending
        if (t.empty()) is.setstate(std::ios::eofbit);
        return is;
      default:
        t += (char) c;
    }
  }
  return is;
}


