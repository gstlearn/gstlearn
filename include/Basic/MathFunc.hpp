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

#include "Matrix/MatrixDense.hpp"

namespace gstlrn
{
class Cheb_Elem;

GSTLEARN_EXPORT Id mvndst_infin(double low, double sup);
GSTLEARN_EXPORT void mvndst(Id n,
                            double *lower,
                            double *upper,
                            Id *infin,
                            double *correl,
                            Id maxpts,
                            double abseps,
                            double releps,
                            double *error,
                            double *value,
                            Id *inform);
GSTLEARN_EXPORT void mvndst2n(const double *lower,
                              const double *upper,
                              const double *means,
                              double *correl,
                              Id maxpts,
                              double abseps,
                              double releps,
                              double *error,
                              double *value,
                              Id *inform);
GSTLEARN_EXPORT void mvndst4(double *lower,
                             double *upper,
                             const double *correl,
                             Id maxpts,
                             double abseps,
                             double releps,
                             double *error,
                             double *value,
                             Id *inform);
GSTLEARN_EXPORT Id besselj_table(double x, double alpha, Id nb, double *b);
GSTLEARN_EXPORT double besselj(double x, Id n);
GSTLEARN_EXPORT Id besselk(double x, double alpha, Id nb, double *bk);
GSTLEARN_EXPORT double loggamma(double parameter);

GSTLEARN_EXPORT double ut_legendre(Id n, double v, bool flagNorm = true);
GSTLEARN_EXPORT VectorDouble ut_legendreVec(Id n, const VectorDouble& vecin, bool flagNorm);
GSTLEARN_EXPORT MatrixDense ut_legendreMatNorm(Id n, const VectorDouble& v);
GSTLEARN_EXPORT MatrixDense ut_legendreAssociatedMat(Id l,
                                                           const VectorDouble &v,
                                                           bool flagNorm = true);

GSTLEARN_EXPORT double ut_flegendre(Id n, Id k0, double theta, bool flagNorm = true);
GSTLEARN_EXPORT double ut_sphericalHarmonic(Id n, Id k, double theta, double phi);
GSTLEARN_EXPORT VectorDouble ut_sphericalHarmonicVec(Id n,
                                                     Id k,
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
GSTLEARN_EXPORT Id ut_chebychev_count(double (*func)(double,
                                                      double,
                                                      const VectorDouble&),
                                       Cheb_Elem *cheb_elem,
                                       double x,
                                       const VectorDouble& blin);
GSTLEARN_EXPORT Id ut_chebychev_coeffs(double (*func)(double,
                                                       double,
                                                       const VectorDouble&),
                                        Cheb_Elem *cheb_elem,
                                        const VectorDouble& blin);
GSTLEARN_EXPORT void ut_vandercorput(Id n,
                                     Id flag_sym,
                                     Id flag_rot,
                                     Id *ntri_arg,
                                     VectorDouble& coord);
GSTLEARN_EXPORT Id ut_icosphere(Id n,
                                 Id flag_rot,
                                 Id* ntri_arg,
                                 VectorDouble& coord);
GSTLEARN_EXPORT double ut_factorial(Id k);
GSTLEARN_EXPORT void ut_log_factorial(Id nbpoly, double *factor);
GSTLEARN_EXPORT MatrixDense* vanDerCorput(Id n, Id nd);
GSTLEARN_EXPORT MatrixDense fillLegendreMatrix(const VectorDouble &r,
                                                     Id legendreOrder);
GSTLEARN_EXPORT Id solve_P2(double a, double b, double c, VectorDouble& x);
GSTLEARN_EXPORT Id solve_P3(double a, double b, double c, double d, VectorDouble& x);
}