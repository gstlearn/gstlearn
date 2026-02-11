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
#define ENUM_CALC_VARIO ECalcVario, VARIOGRAM,                                           \
                        VARIOGRAM, 0, "Variogram",                                       \
                        COVARIANCE, 1, "Covariance",                                     \
                        COVARIOGRAM, 2, "Transitive Covariogram",                        \
                        MADOGRAM, 3, "Madogram",                                         \
                        RODOGRAM, 4, "Rodogram",                                         \
                        POISSON, 5, "Poisson Variogram",                                 \
                        GENERAL1, 6, "Generalized Variogram of order 1",                 \
                        GENERAL2, 7, "Generalized Variogram of order 2",                 \
                        GENERAL3, 8, "Generalized Variogram of order 3",                 \
                        COVARIANCE_NC, 9, "Non-centered Covariance",                     \
                        ORDER4, 10, "Order-4 Variogram",                                 \
                        TRANS1, 11, "Cross-to-Simple Variogram ratio G12/G1",            \
                        TRANS2, 12, "Cross-to-Simple Variogram ratio G12/G2",            \
                        BINORMAL, 13, "Cross-to-Simple Variogram ratio G12/sqrt(G1*G2)", \
                        CORRELOGRAM, 14, "Correlation Variogram"

ENUM_DECLARE(ENUM_CALC_VARIO)

// The next Global variable is there to store additional attributes
// It must be dimensioned parallel to the above ENUM structure

struct qualifier
{
  bool isSymmetric; // True if the tool is symmetric (e.g. Variogram), false otherwise (e.g. Covariance)
  bool isCentered;  // True if this tool is centered by Mean (e.g. Covariance), false otherwise (e.g. Variogram)
  bool isScaled;    // True if this tool is scaled by St. Dev. (e.g. Correlation), false otherwise (e.g. Variogram)
  bool isFittable;  // True if this tool can be used for fitting a model
};

const std::map<std::string, qualifier> ECalcVarioAttr =
  {
    {"UNDEFINED", {true, false, false, true}},
    {"VARIOGRAM", {true, false, false, true}},
    {"COVARIANCE", {false, true, false, true}},
    {"COVARIOGRAM", {false, true, false, true}},
    {"MADOGRAM", {true, false, false, false}},
    {"RODOGRAM", {true, false, false, false}},
    {"POISSON", {true, false, false, true}},
    {"GENERAL1", {true, false, false, false}},
    {"GENERAL2", {true, false, false, false}},
    {"GENERAL3", {true, false, false, false}},
    {"COVARIANCE_NC", {false, true, false, true}},
    {"ORDER4", {true, false, false, true}},
    {"TRANS1", {true, false, false, false}},
    {"TRANS2", {true, false, false, false}},
    {"BINORMAL", {true, false, false, false}},
    {"CORRELOGRAM", {false, true, true, true}}};