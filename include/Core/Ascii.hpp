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

namespace gstlrn
{
GSTLEARN_EXPORT void ascii_study_define(const char* study);
GSTLEARN_EXPORT void ascii_environ_read(String& filename, bool verbose);
GSTLEARN_EXPORT void ascii_filename(const char* type, Id rank, Id mode, String& filename);
GSTLEARN_EXPORT void ascii_simu_read(String& filename, bool verbose, Id* nbsimu, Id* nbtuba, Id* seed);
GSTLEARN_EXPORT bool ascii_option_defined(const String& filename,
                                          const char* option_name,
                                          Id* answer,
                                          bool verbose = false);
} // namespace gstlrn
