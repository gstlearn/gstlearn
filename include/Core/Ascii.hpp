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

GSTLEARN_EXPORT void ascii_study_define(const char* study);
GSTLEARN_EXPORT void ascii_environ_read(char* file_name, int verbose);
GSTLEARN_EXPORT void ascii_filename(const char* type, int rank, int mode, char* filename);
GSTLEARN_EXPORT void ascii_simu_read(char* file_name, int verbose, int* nbsimu, int* nbtuba, int* seed);
GSTLEARN_EXPORT int ascii_option_defined(const char* file_name,
                                         int verbose,
                                         const char* option_name,
                                         int type,
                                         void* answer);
