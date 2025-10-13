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

#include "Basic/VectorNumT.hpp"
#include "Enum/EOperator.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <cmath>
#include <map>

namespace gstlrn
{
typedef struct
{
  Id number;
  Id nvalid;
  double mini;
  double maxi;
  double delta;
  double mean;
  double stdv;
} StatResults;

GSTLEARN_EXPORT bool isInteger(double value, double eps = EPSILON10);
GSTLEARN_EXPORT Id getClosestInteger(double value);
GSTLEARN_EXPORT bool isMultiple(Id nbig, Id nsmall);
GSTLEARN_EXPORT bool isOdd(Id number);
GSTLEARN_EXPORT bool isEven(Id number);
GSTLEARN_EXPORT bool isZero(double value, double eps = EPSILON10);
GSTLEARN_EXPORT bool isOne(double value, double eps = EPSILON10);
GSTLEARN_EXPORT bool isEqual(double v1, double v2, double eps = EPSILON10);
GSTLEARN_EXPORT double getMin(double val1, double val2);
GSTLEARN_EXPORT double getMax(double val1, double val2);
GSTLEARN_EXPORT double ut_deg2rad(double angle);
GSTLEARN_EXPORT double ut_rad2deg(double angle);

GSTLEARN_EXPORT bool isEqualExtended(double v1,
                                     double v2,
                                     double eps           = EPSILON10,
                                     bool flagRelative    = true,
                                     bool flagAbsolute    = false,
                                     const String& string = "");

// No need this stuff through SWIG (because we use target language NAs)
#ifndef SWIG

GSTLEARN_EXPORT bool FFFF(double value); // TODO isNA<double>
GSTLEARN_EXPORT bool IFFFF(Id value);   // TODO isNA<Id>
GSTLEARN_EXPORT double getTEST();        // TODO getNA<double>
GSTLEARN_EXPORT Id getITEST();          // TODO getNA<Id>

#  define DOUBLE_NA TEST
#  define INT_NA    ITEST
#  define STRING_NA "NA"
#  define FLOAT_NA  static_cast<float>(TEST) // 1.234e30 is ok for 4 bytes but needs a cast for Windows

template<typename T>
inline T getNA();
template<>
inline double getNA()
{
  return DOUBLE_NA;
}
template<>
inline Id getNA()
{
  return INT_NA;
}
template<>
inline String getNA()
{
  return STRING_NA;
}
template<>
inline float getNA()
{
  return FLOAT_NA;
}

template<typename T>
inline bool isNA(const T& v);
template<>
inline bool isNA(const double& v)
{
  return (v == getNA<double>() || std::isnan(v) || std::isinf(v));
}
template<>
inline bool isNA(const Id& v)
{
  return (v == getNA<Id>());
}
template<>
inline bool isNA(const String& v)
{
  return (v == getNA<String>());
}
template<>
inline bool isNA(const float& v)
{
  return (v == getNA<float>() || std::isnan(v) || std::isinf(v));
}

#endif // SWIG

// Other Utility functions

GSTLEARN_EXPORT void ut_sort_double(Id safe, Id nech, Id* ind, double* value);
GSTLEARN_EXPORT StatResults ut_statistics(Id nech,
                                          const double* tab,
                                          const double* sel = nullptr,
                                          const double* wgt = nullptr);
GSTLEARN_EXPORT void ut_stats_mima_print(const char* title, Id nech, double* tab, double* sel);
GSTLEARN_EXPORT void ut_facies_statistics(Id nech,
                                          double* tab,
                                          double* sel,
                                          Id* nval,
                                          Id* mini,
                                          Id* maxi);
GSTLEARN_EXPORT void ut_classify(Id nech,
                                 const double* tab,
                                 double* sel,
                                 Id nclass,
                                 double start,
                                 double pas,
                                 Id* nmask,
                                 Id* ntest,
                                 Id* nout,
                                 Id* classe);
GSTLEARN_EXPORT double ut_median(VectorDouble& tab, Id ntab);
GSTLEARN_EXPORT double ut_cnp(Id n, Id k);
GSTLEARN_EXPORT MatrixSquare ut_pascal(Id ndim);
GSTLEARN_EXPORT VectorInt ut_combinations(Id n, Id maxk, Id* ncomb);
GSTLEARN_EXPORT void ut_shuffle_array(Id nrow, Id ncol, VectorDouble& tab);

GSTLEARN_EXPORT VectorInt getListActiveToAbsolute(const VectorDouble& sel);
GSTLEARN_EXPORT std::map<Id, Id> getMapAbsoluteToRelative(const VectorDouble& sel,
                                                            bool verbose = false);
GSTLEARN_EXPORT Id getRankMapAbsoluteToRelative(const std::map<Id, Id>& map, Id iabs);
GSTLEARN_EXPORT Id getRankMapRelativeToAbsolute(const std::map<Id, Id>& map, Id irel);

typedef double (*operate_function)(double);
GSTLEARN_EXPORT operate_function operate_Identify(Id oper);
GSTLEARN_EXPORT double operate_Identity(double x);
GSTLEARN_EXPORT double operate_Inverse(double x);
GSTLEARN_EXPORT double operate_Square(double x);
GSTLEARN_EXPORT double operate_InverseSquare(double x);
GSTLEARN_EXPORT double operate_Sqrt(double x);
GSTLEARN_EXPORT double operate_InverseSqrt(double x);
GSTLEARN_EXPORT double modifyOperator(const EOperator& oper, double oldval, double value);

GSTLEARN_EXPORT double roundZero(double value, double eps = EPSILON6);

GSTLEARN_EXPORT double truncateDecimals(double value, Id ndec = 0);
GSTLEARN_EXPORT double truncateDigits(double value, Id ndigits);

GSTLEARN_EXPORT void print_range(const char* title, Id ntab, const double* tab, const double* sel);

GSTLEARN_EXPORT void convertIndptrToIndices(Id ncumul, const I32* cumul, I32* tab);
} // namespace gstlrn
