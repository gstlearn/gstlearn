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
#

/**
 * TODO : Documentation
 */
#define ENUM_CALC_VARIO ECalcVario, UNDEFINED,                                \
                        UNDEFINED, -1, "Undefined",                           \
                        VARIOGRAM, 0, "Variogram",                            \
                        COVARIANCE, 1, "Covariance",                          \
                        COVARIOGRAM, 2, "Transitive Covariogram",             \
                        MADOGRAM, 3, "Madogram",                              \
                        RODOGRAM, 4, "Rodogram",                              \
                        POISSON, 5, "Poisson Variogram",                      \
                        GENERAL1, 6, "Generalized Variogram of order 1",      \
                        GENERAL2, 7, "Generalized Variogram of order 2",      \
                        GENERAL3, 8, "Generalized Variogram of order 3",      \
                        COVARIANCE_NC, 9, "Non-centered Covariance",          \
                        ORDER4, 10, "Order-4 Variogram",                      \
                        TRANS1, 11, "Cross-to-Simple Variogram ratio G12/G1", \
                        TRANS2, 12, "Cross-to-Simple Variogram ratio G12/G2", \
                        BINORMAL, 13, "Cross-to-Simple Variogram ratio G12/sqrt(G1*G2)"

ENUM_DECLARE(ENUM_CALC_VARIO)

// The next Global variable is there to store additional attributes
// It must be dimensioned parallel to the above ENUM structure

struct qualifier
{
  bool isSymmetric; // True if the tool is symmetric (e.g. Variogram), false otherwise (e.g. Covariance)
  bool isFittable;  // True if this tool can be used for fitting a model
};

inline const std::map<std::string, qualifier> ECalcVarioAttr =
  {
    {"UNDEFINED", {true, true}},
    {"VARIOGRAM", {true, true}},
    {"COVARIANCE", {false, true}},
    {"COVARIOGRAM", {false, true}},
    {"MADOGRAM", {true, false}},
    {"RODOGRAM", {true, false}},
    {"POISSON", {true, true}},
    {"GENERAL1", {false, false}},
    {"GENERAL2", {false, false}},
    {"GENERAL3", {false, false}},
    {"COVARIANCE_NC", {false, true}},
    {"ORDER4", {true, false}},
    {"TRANS1", {false, false}},
    {"TRANS2", {false, false}},
    {"BINORMAL", {false, false}}};
