/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
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
