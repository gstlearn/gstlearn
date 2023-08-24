/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

/**
 * TODO : Documentation
 */
#define ENUM_CALC_VARIO ECalcVario, UNDEFINED,\
                        UNDEFINED,    -1, "Undefined",\
                        VARIOGRAM,     0, "Variogram",\
                        COVARIANCE,    1, "Covariance",\
                        COVARIOGRAM,   2, "Transitive Covariogram",\
                        MADOGRAM,      3, "Madogram",\
                        RODOGRAM,      4, "Rodogram",\
                        POISSON,       5, "Poisson",\
                        GENERAL1,      6, "Generalized Variogram of order 1",\
                        GENERAL2,      7, "Generalized Variogram of order 2",\
                        GENERAL3,      8, "Generalized Variogram of order 3",\
                        COVARIANCE_NC, 9, "Non-centered Covariance",\
                        ORDER4,       10, "Order-4 Variogram",\
                        TRANS1,       11, "Transition probability G12/G1",\
                        TRANS2,       12, "Transition probability G12/G2",\
                        BINORMAL,     13, "Binormal hypothesis G12/sqrt(G1 * G2)"

ENUM_DECLARE(ENUM_CALC_VARIO)
