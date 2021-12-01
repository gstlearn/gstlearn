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

#include "geoslib_define.h"

/// TODO to be removed when VectorDouble will be OK
VectorInt ut_vector_int(int nval, int value = 0.);
VectorDouble ut_vector_double(int nval, double value = 0.);
VectorVectorDouble ut_vector_vector_double(int nval1, int nval2, double value = 0.);
VectorVectorInt ut_vector_vector_int(int nval1, int nval2, int value = 0);
void ut_vector_display(const String& title, const VectorDouble& vect);
void ut_vector_display(const String& title, const VectorVectorDouble& vect);
void ut_vector_display_stats(const String& title, const VectorDouble& vect);
void ut_vector_display_range(const String& title, const VectorDouble& vect);
void ut_ivector_display(const String& title, const VectorInt& vect);
String ut_vector_string(const VectorDouble& vec);
String ut_vector_string(const VectorVectorDouble& vec);
String ut_ivector_string(const VectorInt& vec);
double ut_vector_max(const VectorDouble& vec, bool flagAbs = false);
double ut_vector_min(const VectorDouble& vec, bool flagAbs = false);
double ut_vector_mean(const VectorDouble& vec);
double ut_vector_var(const VectorDouble& vec);
double ut_vector_stdv(const VectorDouble& vec);
double ut_vector_norm(const VectorDouble& vec);
double ut_vector_inner_product(const VectorDouble& vec1,
                               const VectorDouble& vec2);
VectorDouble ut_vector_cross_product(const VectorDouble& vec1,
                                     const VectorDouble& vec2);
bool ut_vector_same(const VectorDouble& v1,
                    const VectorDouble& v2,
                    double eps = EPSILON10);
bool ut_ivector_same(const VectorInt& v1, const VectorInt& v2);
void ut_vector_fill(VectorDouble& vec, double v, int size = 0);
void ut_ivector_fill(VectorInt& vec, int v, int size = 0);
VectorDouble ut_vector_add(const VectorDouble& vec1, const VectorDouble& vec2);
void ut_vector_add_inplace(VectorDouble& vec1, const VectorDouble& vec2);
VectorDouble ut_vector_subtract(const VectorDouble& vec1,
                                const VectorDouble& vec2);
VectorDouble ut_vector_power(const VectorDouble&vec, double power);
void ut_vector_cumul(VectorDouble& vec1, const VectorDouble& vec2, double coeff);
void ut_vector_copy(VectorDouble& vec1, const VectorDouble& vec2);
void ut_vector_multiply_inplace(VectorDouble& vec, double v);
void ut_vector_divide_inplace(VectorDouble& vec, double v);
void ut_vector_addval(VectorDouble& vec, double v);
void ut_vector_sum(const VectorDouble& vec1,const VectorDouble& vec2,VectorDouble& res );
void ut_ivector_addval(VectorInt& vec, int v);
void ut_vector_divide_vec(VectorDouble& vec, const VectorDouble& v);
int  ut_vector_count_undefined(const VectorDouble& vec);

VectorInt    ut_vector_sample(int ntotal, double proportion, int seed = 242141);

VectorDouble ut_vector_simulate_uniform(int n, double mini = 0., double maxi = 1.);
VectorDouble ut_vector_simulate_gaussian(int n, double mean = 0., double sigma = 1.);

int ut_ivector_prod(const VectorInt nx);
VectorInt ut_ivector_sequence(int number, int ideb = 0);
VectorDouble ut_vector_sequence(double valFrom, double valTo, double valStep);

int ut_vector_size(const VectorInt& vec);
int ut_vector_size(const VectorDouble& vec);
int ut_vector_size(const VectorVectorInt& vec);
int ut_vector_size(const VectorVectorDouble& vec);
