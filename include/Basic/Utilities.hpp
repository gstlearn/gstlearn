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

#define DOUBLE_NA TEST
#define INT_NA    ITEST

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

template <typename T> inline T getNAValue()     { return ValueNA<T>::getNA(); }
template <typename T> inline T isNA(const T& v) { return (v == ValueNA<T>::getNA()); }
GSTLEARN_EXPORT bool   isInteger(double value, double eps = EPSILON10);
GSTLEARN_EXPORT int    getClosestInteger(double value);
GSTLEARN_EXPORT bool   isMultiple(int nbig, int nsmall);
GSTLEARN_EXPORT bool   isOdd(int number);
GSTLEARN_EXPORT bool   isEven(int number);
GSTLEARN_EXPORT int    FFFF(double value);
GSTLEARN_EXPORT int    IFFFF(int value);
GSTLEARN_EXPORT double getTEST();
GSTLEARN_EXPORT int    getITEST();
GSTLEARN_EXPORT double getMin(double val1, double val2);
GSTLEARN_EXPORT double getMax(double val1, double val2);

