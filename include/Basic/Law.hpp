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

#include "gstlearn_export.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Basic/VectorNumT.hpp"

GSTLEARN_EXPORT void law_set_old_style(bool style);
GSTLEARN_EXPORT int law_get_random_seed(void);
GSTLEARN_EXPORT void law_set_random_seed(int seed);
GSTLEARN_EXPORT double law_uniform(double mini=0., double maxi=1.);
GSTLEARN_EXPORT int law_int_uniform(int mini, int maxi);
GSTLEARN_EXPORT double law_gaussian(double mean = 0., double sigma = 1.);
GSTLEARN_EXPORT double law_exponential(double lambda = 1.);
GSTLEARN_EXPORT double law_gamma(double alpha, double beta = 1.);
GSTLEARN_EXPORT double law_df_poisson(int i, double parameter);
GSTLEARN_EXPORT VectorDouble law_df_poisson_vec(VectorInt is, double parameter);
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
GSTLEARN_EXPORT double law_df_bigaussian(VectorDouble& vect,
                                         VectorDouble& mean,
                                         MatrixSquareSymmetric& correl);
GSTLEARN_EXPORT double law_df_quadgaussian(VectorDouble &vect,
                                           MatrixSquareSymmetric &correl);
GSTLEARN_EXPORT double law_df_multigaussian(VectorDouble &vect,
                                            MatrixSquareSymmetric& correl);
GSTLEARN_EXPORT VectorInt law_random_path(int nech);
GSTLEARN_EXPORT double* law_exp_sample(const double *tabin,
                                       int mode,
                                       int nvar,
                                       int nechin,
                                       int nechout,
                                       int niter,
                                       int nconst,
                                       double *consts,
                                       int seed,
                                       double percent);
GSTLEARN_EXPORT int sampleInteger(int minit, int maxi);
