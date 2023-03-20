/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
GSTLEARN_EXPORT double ut_legendre(int flag_norm, int n, double v);
GSTLEARN_EXPORT double ut_flegendre(int flag_norm, int n, int k0, double theta);
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
                                                      int,
                                                      double*),
                                       Cheb_Elem *cheb_elem,
                                       double x,
                                       int nblin,
                                       double *blin);
GSTLEARN_EXPORT int ut_chebychev_coeffs(double (*func)(double,
                                                       double,
                                                       int,
                                                       double*),
                                        Cheb_Elem *cheb_elem,
                                        int nblin,
                                        double *blin);
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
