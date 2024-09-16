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

#include "Matrix/MatrixRectangular.hpp"

class Cheb_Elem;

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
GSTLEARN_EXPORT void mvndst2n(const double *lower,
                              const double *upper,
                              const double *means,
                              double *correl,
                              int maxpts,
                              double abseps,
                              double releps,
                              double *error,
                              double *value,
                              int *inform);
GSTLEARN_EXPORT void mvndst4(double *lower,
                             double *upper,
                             const double *correl,
                             int maxpts,
                             double abseps,
                             double releps,
                             double *error,
                             double *value,
                             int *inform);
GSTLEARN_EXPORT int besselj_table(double x, double alpha, int nb, double *b);
GSTLEARN_EXPORT double besselj(double x, int n);
GSTLEARN_EXPORT int besselk(double x, double alpha, int nb, double *bk);
GSTLEARN_EXPORT double loggamma(double parameter);

GSTLEARN_EXPORT double ut_legendre(int n, double v, bool flagNorm = true);
GSTLEARN_EXPORT VectorDouble ut_legendreVec(int n, const VectorDouble& vecin, bool flagNorm);
GSTLEARN_EXPORT MatrixRectangular ut_legendreMatNorm(int n, const VectorDouble& v);
GSTLEARN_EXPORT MatrixRectangular ut_legendreAssociatedMat(int l,
                                                           const VectorDouble &v,
                                                           bool flagNorm = true);

GSTLEARN_EXPORT double ut_flegendre(int n, int k0, double theta, bool flagNorm = true);
GSTLEARN_EXPORT double ut_sphericalHarmonic(int n, int k, double theta, double phi);
GSTLEARN_EXPORT VectorDouble ut_sphericalHarmonicVec(int n,
                                                     int k,
                                                     VectorDouble theta,
                                                     VectorDouble phi);
GSTLEARN_EXPORT double golden_search(double (*func_evaluate)(double test,
                                                             void *user_data),
                                                             void *user_data,
                                                             double tolstop,
                                                             double a0,
                                                             double c0,
                                                             double *test_loc,
                                                             double *niter);
GSTLEARN_EXPORT int ut_chebychev_count(double (*func)(double,
                                                      double,
                                                      const VectorDouble&),
                                       Cheb_Elem *cheb_elem,
                                       double x,
                                       const VectorDouble& blin);
GSTLEARN_EXPORT int ut_chebychev_coeffs(double (*func)(double,
                                                       double,
                                                       const VectorDouble&),
                                        Cheb_Elem *cheb_elem,
                                        const VectorDouble& blin);
GSTLEARN_EXPORT void ut_vandercorput(int n,
                                     int flag_sym,
                                     int flag_rot,
                                     int *ntri_arg,
                                     double **coor_arg);
GSTLEARN_EXPORT int ut_icosphere(int n, int flag_rot, int *ntri_arg, double **coor_arg);
GSTLEARN_EXPORT double ut_factorial(int k);
GSTLEARN_EXPORT void ut_log_factorial(int nbpoly, double *factor);
GSTLEARN_EXPORT MatrixRectangular* vanDerCorput(int n, int nd);
GSTLEARN_EXPORT MatrixRectangular fillLegendreMatrix(const VectorDouble &r,
                                                     int legendreOrder);
GSTLEARN_EXPORT int solve_P2(double a, double b, double c, VectorDouble& x);
GSTLEARN_EXPORT int solve_P3(double a, double b, double c, double d, VectorDouble& x);
