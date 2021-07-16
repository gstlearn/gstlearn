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

int mvndst_infin(double low, double sup);
void mvndst(int n,
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
void mvndst2n(double *lower,
              double *upper,
              double *means,
              double *correl,
              int maxpts,
              double abseps,
              double releps,
              double *error,
              double *value,
              int *inform);
void mvndst4(double *lower,
             double *upper,
             double *correl,
             int maxpts,
             double abseps,
             double releps,
             double *error,
             double *value,
             int *inform);
int bessel_j(double x, double alpha, int nb, double *b);
int bessel_k(double x, double alpha, int nb, double *bk);
double loggamma(double parameter);
