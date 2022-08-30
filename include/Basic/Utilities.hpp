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
#include "Geometry/Geometry.hpp"

GSTLEARN_EXPORT bool   isInteger(double value, double eps = EPSILON10);
GSTLEARN_EXPORT int    getClosestInteger(double value);
GSTLEARN_EXPORT bool   isMultiple(int nbig, int nsmall);
GSTLEARN_EXPORT bool   isOdd(int number);
GSTLEARN_EXPORT bool   isEven(int number);
GSTLEARN_EXPORT int    FFFF(double value); // TODO isNA<double>
GSTLEARN_EXPORT int    IFFFF(int value);   // TODO isNA<int.
GSTLEARN_EXPORT double getTEST();  // TODO getNAValue<double>
GSTLEARN_EXPORT int    getITEST(); // TODO getNAValue<int>
GSTLEARN_EXPORT double getMin(double val1, double val2);
GSTLEARN_EXPORT double getMax(double val1, double val2);
GSTLEARN_EXPORT double ut_deg2rad(double angle);
GSTLEARN_EXPORT double ut_rad2deg(double angle);

#define DOUBLE_NA TEST
#define INT_NA    ITEST
#define STRING_NA "NA"    // TODO search for this string and replace

template<typename T> class ValueNA;

// Define NA value for double
template <> class ValueNA<double>
{
public:
  static inline double getNA() { return DOUBLE_NA; }
};

// Define NA value for int
template <> class ValueNA<int>
{
public:
    static inline int getNA() { return INT_NA; }
};

// Define NA value for String
template <> class ValueNA<String>
{
public:
  static inline String getNA() { return STRING_NA; }
};

template <typename T> inline T    getNAValue()     { return ValueNA<T>::getNA(); }
template <typename T> inline bool isNA(const T& v) { return (v == ValueNA<T>::getNA()); }

// Other Utiity functions


GSTLEARN_EXPORT void ut_tab_unique(int ntab, double *tab, int *neff);
GSTLEARN_EXPORT void ut_sort_double(int safe, int nech, int *ind, double *value);
GSTLEARN_EXPORT void ut_sort_int(int safe, int nech, int *ind, int *value);
GSTLEARN_EXPORT void ut_statistics(int nech,
                                   double *tab,
                                   double *sel,
                                   double *wgt,
                                   int *nval,
                                   double *mini,
                                   double *maxi,
                                   double *delta,
                                   double *mean,
                                   double *stdv);
GSTLEARN_EXPORT void ut_stats_mima(int nech,
                   double *tab,
                   double *sel,
                   int *nvalid,
                   double *mini,
                   double *maxi);
GSTLEARN_EXPORT void ut_stats_mima_print(const char *title, int nech, double *tab, double *sel);
GSTLEARN_EXPORT void ut_facies_statistics(int nech,
                                          double *tab,
                                          double *sel,
                                          int *nval,
                                          int *mini,
                                          int *maxi);
GSTLEARN_EXPORT void ut_classify(int nech,
                                 double *tab,
                                 double *sel,
                                 int nclass,
                                 double start,
                                 double pas,
                                 int *nmask,
                                 int *ntest,
                                 int *nout,
                                 int *classe);
GSTLEARN_EXPORT double ut_median(double *tab, int ntab);
GSTLEARN_EXPORT void ut_normalize(int ntab, double *tab);
GSTLEARN_EXPORT double ut_cnp(int n, int k);
GSTLEARN_EXPORT double* ut_pascal(int ndim);
GSTLEARN_EXPORT int* ut_combinations(int n, int maxk, int *ncomb);
GSTLEARN_EXPORT void ut_shuffle_array(int nrow, int ncol, double *tab);

