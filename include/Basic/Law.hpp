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

GSTLEARN_EXPORT void law_set_old_style(bool style);
GSTLEARN_EXPORT int law_get_random_seed(void);
GSTLEARN_EXPORT void law_set_random_seed(int seed);
GSTLEARN_EXPORT double law_uniform(double mini=0., double maxi=1.);
GSTLEARN_EXPORT int law_int_uniform(int mini, int maxi);
GSTLEARN_EXPORT double law_gaussian(double mean = 0., double sigma = 1.);
GSTLEARN_EXPORT double law_exponential(double lambda = 1.);
GSTLEARN_EXPORT double law_gamma(double alpha, double beta = 1.);
GSTLEARN_EXPORT int law_poisson(double parameter);
GSTLEARN_EXPORT double law_stable_standard_agd(double alpha, double beta);
GSTLEARN_EXPORT double law_stable_standard_a1gd(double beta);
GSTLEARN_EXPORT double law_stable_standard_abgd(double alpha);
GSTLEARN_EXPORT double law_stable_a(double alpha,
                                    double beta,
                                    double gamma,
                                    double delta);
GSTLEARN_EXPORT double law_stable_a1(double beta, double gamma, double delta);
GSTLEARN_EXPORT double law_stable(double alpha,
                                  double beta,
                                  double gamma,
                                  double delta);
GSTLEARN_EXPORT int law_binomial(int n, double p);
GSTLEARN_EXPORT double law_beta1(double parameter1, double parameter2);
GSTLEARN_EXPORT double law_beta2(double parameter1, double parameter2);
GSTLEARN_EXPORT double law_df_gaussian(double value);
GSTLEARN_EXPORT double law_dnorm(double value, double mean, double std);
GSTLEARN_EXPORT double law_cdf_gaussian(double value);
GSTLEARN_EXPORT double law_invcdf_gaussian(double value);
GSTLEARN_EXPORT double law_gaussian_between_bounds(double binf, double bsup);
GSTLEARN_EXPORT double law_df_bigaussian(double *vect,
                                         double *mean,
                                         double *corr);
GSTLEARN_EXPORT double law_df_quadgaussian(double *vect, double *corr);
GSTLEARN_EXPORT double law_df_multigaussian(int nvar,
                                            double *vect,
                                            double *corr);
GSTLEARN_EXPORT void law_random_path(int nech, int *path);
GSTLEARN_EXPORT double* law_exp_sample(double *tabin,
                                       int mode,
                                       int nvar,
                                       int nechin,
                                       int nechout,
                                       int niter,
                                       int nconst,
                                       double *consts,
                                       int seed,
                                       double percent);
