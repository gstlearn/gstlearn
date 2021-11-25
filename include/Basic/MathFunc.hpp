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

#include "gstlearn_export.hpp"

GSTLEARN_EXPORT int mvndst_infin(double low, double sup);
GSTLEARN_EXPORT void mvndst(int n,
                            double *lower,
                            double *upper,
                            int *infin,
                            double *correl,
                            int maxpts,
                            double abseps,
                            double releps,
                            double *error,
                            double *value,
                            int *inform);
GSTLEARN_EXPORT void mvndst2n(double *lower,
                              double *upper,
                              double *means,
                              double *correl,
                              int maxpts,
                              double abseps,
                              double releps,
                              double *error,
                              double *value,
                              int *inform);
GSTLEARN_EXPORT void mvndst4(double *lower,
                             double *upper,
                             double *correl,
                             int maxpts,
                             double abseps,
                             double releps,
                             double *error,
                             double *value,
                             int *inform);
GSTLEARN_EXPORT int bessel_j(double x, double alpha, int nb, double *b);
GSTLEARN_EXPORT int bessel_k(double x, double alpha, int nb, double *bk);
GSTLEARN_EXPORT double loggamma(double parameter);
