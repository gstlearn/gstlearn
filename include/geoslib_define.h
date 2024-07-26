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
#include "Basic/WarningMacro.hpp"
#include "Basic/RepeatMacro.hpp"

// WARNING: Make this include list as small as possible!
#include <string>
typedef std::string String;
typedef unsigned char UChar;

#define EPSILON1   1.e-1
#define EPSILON2   1.e-2
#define EPSILON3   1.e-3
#define EPSILON4   1.e-4
#define EPSILON5   1.e-5
#define EPSILON6   1.e-6
#define EPSILON7   1.e-7
#define EPSILON8   1.e-8
#define EPSILON9   1.e-9
#define EPSILON10  1.e-10
#define EPSILON13  1.e-13
#define EPSILON20  1.e-20

#define EPSGRAD    1.e-5

// Macro for preventing warning : unused variable.
// To be used like: DECLARE_UNUSED(a, b, c)
#define DECLARE_UNUSED_(x) (void)x;
#define DECLARE_UNUSED(...) EXPAND(REPEAT(DECLARE_UNUSED_, __VA_ARGS__))

// Declare the function which has a specific implementation
// in the Target language.
// This function must be:
// - declared in rgstlearn.i or pygstlearn.i
// - be called as 'classname'_toTL
#define DECLARE_TOTL void toTL() const {};

// No need to this stuff through SWIG (using target language NAs)
// => Not really : Using customized SWIG 4.2.0b, TEST is often a default argument value!
//#ifndef SWIG
#define TEST      1.234e30
#define TEST_COMP 1.000e30
#define ITEST    -1234567
//#endif

#define ASCII_TEST    -999.

#define BUFFER_LENGTH 10000
#define STRING_LENGTH   100
#define LOCAL_SIZE       10
#define LONG_SIZE     10000
#define GV_PI  3.14159265358979323846264338328
#define GV_EE  2.732
#define MIN(a,b)       (((a) < (b)) ?  (a) : (b))
#define MAX(a,b)       (((a) > (b)) ?  (a) : (b))
#define ABS(a)         (((a) <  0.) ? -(a) : (a))
#define SIGN(s,a)      (((s) <  0.) ? -(a) : (a))
#define M_R(tab,n,i,j) (tab[(n) * (i) + (j)])
#define IS_GAUSS_DEF(x)     (x > THRESH_INF  && x <  THRESH_SUP)
#define ISNOT_GAUSS_DEF(x)  (x <= THRESH_INF || x >= THRESH_SUP)

#define MAX_INT     1000000000
#define MAX_PARAM   1000
#define MAX_EXP     5      // Maximum value for exp(-h)
#define MAX_EXP2    10     // Maximum value for exp(-h^2)

#define THRESH_INF      -10
#define THRESH_SUP       10

#define EARTH_RADIUS  6371.

// Hide warnings C4251 under windows: https://stackoverflow.com/a/22054743
#ifndef SWIG
DISABLE_WARNING_NOT_EXPORTED_FROM_DLL
DISABLE_WARNING_BASE_NOT_EXPORTED_FROM_DLL
#endif

