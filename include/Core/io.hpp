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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <cstdarg>
#include <cstdio>

class VectorUChar;

namespace gstlrn
{
GSTLEARN_EXPORT Id _file_read(FILE* file, const char* format, va_list ap);
GSTLEARN_EXPORT Id _file_get_ncol(FILE* file);
GSTLEARN_EXPORT void _token_delimitors(const char del_com, const char del_sep, const char del_blk);
GSTLEARN_EXPORT FILE* _file_open(const char* filename, Id mode);
GSTLEARN_EXPORT Id _record_read(FILE* file, const char* format, void* out);
GSTLEARN_EXPORT void _file_write(FILE* file, const char* format, va_list ap);
GSTLEARN_EXPORT void _buffer_write(String& buffer, const char* format, va_list ap);
GSTLEARN_EXPORT void _lire_string(const char* question,
                                  Id flag_def,
                                  const char* valdef,
                                  char* answer);
GSTLEARN_EXPORT Id _lire_int(const char* question, Id flag_def, Id valdef, Id valmin, Id valmax);
GSTLEARN_EXPORT double _lire_double(const char* question,
                                    Id flag_def,
                                    double valdef,
                                    double valmin,
                                    double valmax);
GSTLEARN_EXPORT Id _lire_logical(const char* question, Id flag_def, Id valdef);
GSTLEARN_EXPORT void _erase_current_string(void);
} // namespace gstlrn
