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

#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixSymmetric.hpp"

namespace gstlrn
{
GSTLEARN_EXPORT void law_set_old_style(bool style);
GSTLEARN_EXPORT Id law_get_random_seed(void);
GSTLEARN_EXPORT void law_set_random_seed(Id seed);
GSTLEARN_EXPORT double law_uniform(double mini = 0., double maxi = 1.);
GSTLEARN_EXPORT Id law_int_uniform(Id mini, Id maxi);
GSTLEARN_EXPORT double law_gaussian(double mean = 0., double sigma = 1.);
GSTLEARN_EXPORT double law_exponential(double lambda = 1.);
GSTLEARN_EXPORT double law_gamma(double alpha, double beta = 1.);
GSTLEARN_EXPORT double law_df_poisson(Id i, double parameter);
GSTLEARN_EXPORT VectorDouble law_df_poisson_vec(VectorInt is, double parameter);
GSTLEARN_EXPORT Id law_poisson(double parameter);
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
GSTLEARN_EXPORT Id law_binomial(Id n, double p);
GSTLEARN_EXPORT double law_beta1(double parameter1, double parameter2);
GSTLEARN_EXPORT double law_beta2(double parameter1, double parameter2);
GSTLEARN_EXPORT double law_df_gaussian(double value);
GSTLEARN_EXPORT double law_dnorm(double value, double mean, double std);
GSTLEARN_EXPORT double law_cdf_gaussian(double value);
GSTLEARN_EXPORT double law_invcdf_gaussian(double value);
GSTLEARN_EXPORT double law_gaussian_between_bounds(double binf, double bsup);
GSTLEARN_EXPORT double law_df_bigaussian(VectorDouble& vect,
                                         VectorDouble& mean,
                                         MatrixSymmetric& correl);
GSTLEARN_EXPORT double law_df_quadgaussian(VectorDouble& vect,
                                           MatrixSymmetric& correl);
GSTLEARN_EXPORT double law_df_multigaussian(VectorDouble& vect,
                                            MatrixSymmetric& correl);
GSTLEARN_EXPORT VectorInt law_random_path(Id nech);
GSTLEARN_EXPORT VectorDouble law_exp_sample(const double* tabin,
                                            Id mode,
                                            Id nvar,
                                            Id nechin,
                                            Id nechout,
                                            Id niter,
                                            Id nconst,
                                            double* consts,
                                            Id seed,
                                            double percent);
GSTLEARN_EXPORT Id sampleInteger(Id minit, Id maxi);
} // namespace gstlrn