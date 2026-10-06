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

#include "Enum/AEnum.hpp"

// clang-format off
#define ENUM_ROLE                                                      \
  ERole, Z,                                                            \
  UNDEFINED, -1, "No locator",                                         \
  X, 0, "Coordinate",                                                  \
  Z, 1, "Variable",                                                    \
  V, 2, "Variance of measurement error",                               \
  F, 3, "External Drift",                                              \
  G, 4, "Gradient component",                                          \
  L, 5, "Lower bound of an inequality",                                \
  U, 6, "Upper bound of an inequality",                                \
  P, 7, "Proportion",                                                  \
  W, 8, "Weight",                                                      \
  C, 9, "Code",                                                        \
  SEL, 10, "Selection",                                                \
  DOM, 11, "Domain",                                                   \
  BLEX, 12, "Block Extension",                                         \
  ADIR, 13, "Dip direction Angle",                                     \
  ADIP, 14, "Dip Angle",                                               \
  SIZE, 15, "Object height",                                           \
  BU, 16, "Fault UP termination",                                      \
  BD, 17, "Fault DOWN termination",                                    \
  TIME, 18, "Time variable",                                           \
  LAYER, 19, "Layer rank",                                             \
  NOSTAT, 20, "Non-stationary parameter",                              \
  TGT, 21, "Tangent",                                                  \
  SIMU, 22, "Conditional or non-conditional simulations",              \
  FACIES, 23, "Facies simulated",                                      \
  GAUSFAC, 24, "Gaussian value for Facies",                            \
  DATE, 25, "Date",                                                    \
  RKLOW, 26, "Rank for lower bound (when discretized)",                \
  RKUP, 27, "Rank for upper bound (when discretized)",                 \
  SUM, 28, "Constraints on the Sum"
// clang-format on

ENUM_DECLARE(ENUM_ROLE)

// The next Global variable is there to store additional attributes
// It must be dimensioned parallel to the above ENUM structure

struct roleQualifier
{
  bool isUnique; // Unicity of this Role
  bool isMultiVersionAvailable; // Multi-version of this Role is available
};

// clang-format off
inline const std::map<std::string_view, roleQualifier> ERoleAttr = {
  {"UNDEFINED", {false, false}},
  {"X", {false, false}},
  {"Z", {false, true}},
  {"V", {false, false}},
  {"F", {false, false}},
  {"G", {false, false}},
  {"L", {false, false}},
  {"U", {false, false}},
  {"P", {false, false}},
  {"W", {true, false}},
  {"C", {true, false}},
  {"SEL", {true, false}},
  {"DOM", {true, false}},
  {"BLEX", {false, false}},
  {"ADIR", {true, false}},
  {"ADIP", {true, false}},
  {"SIZE", {true, false}},
  {"BU", {true, false}},
  {"BD", {true, false}},
  {"TIME", {false, false}},
  {"LAYER", {true, false}},
  {"NOSTAT", {false, false}},
  {"TGT", {false, false}},
  {"SIMU", {false, true}},
  {"FACIES", {false, true}},
  {"GAUSFAC", {false, true}},
  {"DATE", {true, false}},
  {"RKLOW", {false, false}},
  {"RKUP", {false, false}},
  {"SUM", {false, false}}
};
// clang-format on
