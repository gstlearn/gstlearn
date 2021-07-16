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

#include <vector>

int law_get_random_seed(void);
void law_set_random_seed(int seed);
double law_uniform(double mini, double maxi);
int law_int_uniform(int mini, int maxi);
double law_gaussian(void);
double law_exponential(void);
double law_gamma(double parameter);
int law_poisson(double parameter);
double law_stable_standard_agd(double alpha, double beta);
double law_stable_standard_a1gd(double beta);
double law_stable_standard_abgd(double alpha);
double law_stable_a(double alpha, double beta, double gamma, double delta);
double law_stable_a1(double beta, double gamma, double delta);
double law_stable(double alpha, double beta, double gamma, double delta);
int law_binomial(int n, double p);
double law_beta1(double parameter1, double parameter2);
double law_beta2(double parameter1, double parameter2);
double law_df_gaussian(double value);
double law_dnorm(double value, double mean, double std);
double law_cdf_gaussian(double value);
double law_invcdf_gaussian(double value);
double law_gaussian_between_bounds(double binf, double bsup);
double law_df_bigaussian(double *vect, double *mean, double *corr);
double law_df_quadgaussian(double *vect, double *corr);
double law_df_multigaussian(int nvar, double *vect, double *corr);
void law_random_path(int nech, int *path);
double *law_exp_sample(double *tabin,
                       int mode,
                       int nvar,
                       int nechin,
                       int nechout,
                       int niter,
                       int nconst,
                       double *consts,
                       int seed,
                       double percent);
