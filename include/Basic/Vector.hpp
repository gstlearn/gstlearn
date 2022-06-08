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
#include "geoslib_define.h"

/// TODO to be removed when VectorDouble will be OK
GSTLEARN_EXPORT VectorInt ut_vector_int(int nval, int value = 0.);
GSTLEARN_EXPORT VectorDouble ut_vector_double(int nval, double value = 0.);
GSTLEARN_EXPORT VectorVectorDouble ut_vector_vector_double(int nval1, int nval2, double value = 0.);
GSTLEARN_EXPORT VectorVectorInt ut_vector_vector_int(int nval1, int nval2, int value = 0);
GSTLEARN_EXPORT VectorInt ut_ivector_set(int* values, int number);
GSTLEARN_EXPORT VectorDouble ut_vector_set(double* values, int number);
GSTLEARN_EXPORT VectorVectorDouble ut_vector_vector_set(double* value, int n1, int n2);
GSTLEARN_EXPORT void ut_vector_display(const String &title,
                                       const VectorDouble& vect);
GSTLEARN_EXPORT void ut_vector_display(const String &title,
                                       const VectorVectorDouble& vect);
GSTLEARN_EXPORT void ut_vector_display(const String &title,
                                       const VectorString& vect);
GSTLEARN_EXPORT void ut_vector_display_stats(const String &title,
                                             const VectorDouble& vect);
GSTLEARN_EXPORT void ut_vector_display_range(const String &title,
                                             const VectorDouble& vect);
GSTLEARN_EXPORT void ut_ivector_display(const String &title,
                                        const VectorInt& vect);
GSTLEARN_EXPORT String ut_vector_string(const VectorDouble& vec);
GSTLEARN_EXPORT String ut_vector_string(const VectorString& vec);
GSTLEARN_EXPORT String ut_ivector_string(const VectorInt& vec);
GSTLEARN_EXPORT int ut_ivector_max(const VectorInt& vec);
GSTLEARN_EXPORT int ut_ivector_min(const VectorInt& vec);
GSTLEARN_EXPORT double ut_vector_max(const VectorDouble& vec,
                                     bool flagAbs = false);
GSTLEARN_EXPORT double ut_vector_min(const VectorDouble& vec,
                                     bool flagAbs = false);
GSTLEARN_EXPORT double ut_vector_cumul(const VectorDouble& vec);
GSTLEARN_EXPORT double ut_vector_mean(const VectorDouble& vec);
GSTLEARN_EXPORT double ut_vector_var(const VectorDouble& vec);
GSTLEARN_EXPORT double ut_vector_stdv(const VectorDouble& vec);
GSTLEARN_EXPORT double ut_vector_norm(const VectorDouble& vec);
GSTLEARN_EXPORT double ut_vector_inner_product(const VectorDouble& vec1,
                                               const VectorDouble& vec2);
GSTLEARN_EXPORT VectorDouble ut_vector_cross_product(const VectorDouble& vec1,
                                                     const VectorDouble& vec2);
GSTLEARN_EXPORT bool ut_vector_constant(const VectorDouble& vect, double refval = TEST);
GSTLEARN_EXPORT bool ut_ivector_constant(const VectorInt& vect, int refval = ITEST);
GSTLEARN_EXPORT bool ut_vector_same(const VectorDouble& v1,
                                    const VectorDouble& v2,
                                    double eps = EPSILON10);
GSTLEARN_EXPORT bool ut_ivector_same(const VectorInt& v1, const VectorInt& v2);
GSTLEARN_EXPORT void ut_vector_fill(VectorDouble& vec, double v, int size = 0);
GSTLEARN_EXPORT void ut_ivector_fill(VectorInt& vec, int v, int size = 0);
GSTLEARN_EXPORT VectorDouble ut_vector_concatenate(const VectorDouble& vec1,
                                                   const VectorDouble& vec2);
GSTLEARN_EXPORT VectorDouble ut_vector_add(const VectorDouble& vec1,
                                           const VectorDouble& vec2);
GSTLEARN_EXPORT void ut_vector_add_inplace(VectorDouble& dest,
                                           const VectorDouble& src);
GSTLEARN_EXPORT VectorDouble ut_vector_subtract(const VectorDouble& vec1,
                                                const VectorDouble& vec2);
GSTLEARN_EXPORT VectorDouble ut_vector_power(const VectorDouble& vec,
                                             double power);
GSTLEARN_EXPORT void ut_vector_cumul(VectorDouble& vec1,
                                     const VectorDouble& vec2,
                                     double coeff);

GSTLEARN_EXPORT std::pair<double,double> ut_vector_rangeVals(const VectorDouble& vec);
GSTLEARN_EXPORT void ut_vector_copy(VectorDouble& vec1,
                                    const VectorDouble& vec2);
GSTLEARN_EXPORT void ut_vector_multiply_inplace(VectorDouble& vec, double v);
GSTLEARN_EXPORT void ut_vector_divide_inplace(VectorDouble& vec, double v);
GSTLEARN_EXPORT VectorDouble ut_vector_inverse(const VectorDouble& vec);

GSTLEARN_EXPORT void ut_vector_addval(VectorDouble& vec, double v);
GSTLEARN_EXPORT void ut_vector_sum(const VectorDouble& vec1,
                                   const VectorDouble& vec2,
                                   VectorDouble& res);
GSTLEARN_EXPORT void ut_ivector_addval(VectorInt& vec, int v);
GSTLEARN_EXPORT void ut_vector_divide_vec(VectorDouble& vec,
                                          const VectorDouble& v);
GSTLEARN_EXPORT int ut_vector_count_undefined(const VectorDouble& vec);

GSTLEARN_EXPORT VectorInt ut_vector_sample(int ntotal,
                                           double proportion = 0.,
                                           int number = 0,
                                           int seed = 242141);

GSTLEARN_EXPORT VectorDouble ut_vector_simulate_uniform(int n,
                                                        double mini = 0.,
                                                        double maxi = 1.);
GSTLEARN_EXPORT VectorDouble ut_vector_simulate_gaussian(int n,
                                                         double mean = 0.,
                                                         double sigma = 1.);
GSTLEARN_EXPORT VectorDouble ut_vector_simulate_bernoulli(int n,
                                                          double proba,
                                                          double vone = 1.,
                                                          double velse = 0.);

GSTLEARN_EXPORT void ut_vector_simulate_gaussian_inplace(VectorDouble & vect,
                                                         double mean = 0.,
                                                         double sigma = 1.);
GSTLEARN_EXPORT int ut_vector_prod(const VectorInt& nx);
GSTLEARN_EXPORT double ut_vector_prod(const VectorDouble& nx);
GSTLEARN_EXPORT VectorInt ut_ivector_sequence(int number, int ideb = 0);
GSTLEARN_EXPORT VectorDouble ut_vector_sequence(double valFrom,
                                                double valTo,
                                                double valStep);

GSTLEARN_EXPORT int ut_vector_size(const VectorInt& vec);
GSTLEARN_EXPORT int ut_vector_size(const VectorDouble& vec);
GSTLEARN_EXPORT int ut_vector_size(const VectorVectorInt& vec);
GSTLEARN_EXPORT int ut_vector_size(const VectorVectorDouble& vec);

GSTLEARN_EXPORT VectorInt ut_ivector_sort(const VectorInt& vecin, bool ascending = true);
GSTLEARN_EXPORT VectorDouble ut_vector_sort(const VectorDouble& vecin, bool ascending = true);
GSTLEARN_EXPORT VectorInt ut_vector_sort_indices(const VectorDouble& vecin);

GSTLEARN_EXPORT double ut_vector_extension_diagonal(const VectorDouble& mini,
                                                    const VectorDouble& maxi);
GSTLEARN_EXPORT VectorDouble ut_vector_angle_from_codir(const VectorDouble& codir);
