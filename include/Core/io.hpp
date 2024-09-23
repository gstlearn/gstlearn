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

#include <stdarg.h>
#include <stdio.h>

GSTLEARN_EXPORT int _file_read(FILE* file, const char* format, va_list ap);
GSTLEARN_EXPORT int _file_get_ncol(FILE* file);
GSTLEARN_EXPORT void _file_delimitors(char del_com, const char* del_sep, char del_blk);
GSTLEARN_EXPORT FILE* _file_open(const char* filename, int mode);
GSTLEARN_EXPORT int _record_read(FILE* file, const char* format, ...);
GSTLEARN_EXPORT int _buffer_read(char** buffer, const char* format, va_list ap);
GSTLEARN_EXPORT void _file_write(FILE* file, const char* format, va_list ap);
GSTLEARN_EXPORT void _buffer_write(char* buffer, const char* format, va_list ap);
GSTLEARN_EXPORT void _lire_string(const char* question,
                                  int flag_def,
                                  const char* valdef,
                                  char* answer);
GSTLEARN_EXPORT int _lire_int(const char* question, int flag_def, int valdef, int valmin, int valmax);
GSTLEARN_EXPORT double _lire_double(const char* question,
                                    int flag_def,
                                    double valdef,
                                    double valmin,
                                    double valmax);
GSTLEARN_EXPORT int _lire_logical(const char* question, int flag_def, int valdef);
GSTLEARN_EXPORT void _erase_current_string(void);