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

#include <Basic/VectorNumT.hpp>

namespace gstlrn
{

GSTLEARN_EXPORT void set_keypair(
  const char* keyword, Id origin, Id nrow, Id ncol, const double* values);
GSTLEARN_EXPORT void app_keypair(
  const char* keyword, Id origin, Id nrow, Id ncol, double* values);
GSTLEARN_EXPORT void set_keypair_int(
  const char* keyword, Id origin, Id nrow, Id ncol, Id* values);
GSTLEARN_EXPORT void app_keypair_int(
  const char* keyword, Id origin, Id nrow, Id ncol, Id* values);
GSTLEARN_EXPORT double get_keypone(const char* keyword, double valdef);
GSTLEARN_EXPORT Id
get_keypair(const char* keyword, Id* nrow, Id* ncol, VectorDouble& values);
GSTLEARN_EXPORT Id
get_keypair_int(const char* keyword, Id* nrow, Id* ncol, VectorInt& values);
GSTLEARN_EXPORT void del_keypair(const char* keyword, Id flag_exact);
GSTLEARN_EXPORT void print_keypair(Id flag_short);
}
