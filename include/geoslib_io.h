/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/******************************************************************************/
#pragma once

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
GSTLEARN_EXPORT void record_close(void);
#ifndef SWIG
GSTLEARN_EXPORT void redefine_message(void (*write_func)(const char*));
GSTLEARN_EXPORT void redefine_error(void (*warn_func)(const char*));
GSTLEARN_EXPORT void redefine_read(void (*read_func)(const char*, char*));
GSTLEARN_EXPORT void redefine_exit(void (*exit_func)(void));
#endif
GSTLEARN_EXPORT void mem_error(Id nbyte);

GSTLEARN_EXPORT void message_extern(const char* string);
GSTLEARN_EXPORT void exit_extern();

GSTLEARN_EXPORT void string_strip_blanks(char* string, Id flag_lead);
GSTLEARN_EXPORT void string_strip_quotes(char* string);

#if defined(_WIN32) || defined(_WIN64)
GSTLEARN_EXPORT char* strsep(char** stringp, const char* delim);
#endif
GSTLEARN_EXPORT void print_current_line(void);

} // namespace gstlrn