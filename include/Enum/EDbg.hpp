/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_DEBUG EDbg, DB, \
                     INTERFACE,   0, "Communication with interface", \
                     DB,          1, "Data Base Management", \
                     NBGH,        2, "Neighborhood Management", \
                     MODEL,       3, "Model Management", \
                     KRIGING,     4, "Kriging Operations", \
                     SIMULATE,    5, "Simulations", \
                     RESULTS,     6, "Kriging Results", \
                     VARIOGRAM,   7, "Variogram calculations", \
                     CONVERGE,    8, "Convergence test", \
                     CONDEXP,     9, "Conditional Expectation", \
                     BAYES,      10, "Bayesian Estimation", \
                     MORPHO,     11, "Morphological operations", \
                     PROPS,      12, "Proportions or Intensities", \
                     UPSCALE,    13, "Upscaling", \
                     SPDE,       14, "S.P.D.E. calculations"

ENUM_DECLARE(ENUM_DEBUG)
