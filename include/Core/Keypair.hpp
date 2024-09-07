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

GSTLEARN_EXPORT void set_keypair(
  const char* keyword, int origin, int nrow, int ncol, const double* values);
GSTLEARN_EXPORT void app_keypair(
  const char* keyword, int origin, int nrow, int ncol, double* values);
GSTLEARN_EXPORT void set_keypair_int(
  const char* keyword, int origin, int nrow, int ncol, int* values);
GSTLEARN_EXPORT void app_keypair_int(
  const char* keyword, int origin, int nrow, int ncol, int* values);
GSTLEARN_EXPORT double get_keypone(const char* keyword, double valdef);
GSTLEARN_EXPORT int
get_keypair(const char* keyword, int* nrow, int* ncol, double** values);
GSTLEARN_EXPORT int
get_keypair_int(const char* keyword, int* nrow, int* ncol, int** values);
GSTLEARN_EXPORT void del_keypair(const char* keyword, int flag_exact);
GSTLEARN_EXPORT void print_keypair(int flag_short);
