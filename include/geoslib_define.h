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

#include "Basic/WarningMacro.hpp"

// WARNING: Make this include list as small as possible!
#include <vector>
#include <string>

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
#define EPSILON20  1.e-20

#define EPSGRAD    1.e-5

#define TEST      1.234e30
#define TEST_COMP 1.000e30
#define ITEST    -1234567

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
#define IS_GAUSS_DEF(x) (x > THRESH_INF && x < THRESH_SUP)

#define MAX_INT     1000000000
#define MAX_PARAM   1000
#define MAX_EXP     5      // Maximum value for exp(-h)
#define MAX_EXP2    10     // Maximum value for exp(-h^2)

#define THRESH_INF      -10
#define THRESH_SUP       10

// Hide warnings C4251 under windows: https://stackoverflow.com/a/22054743
#ifndef SWIG
DISABLE_WARNING_NOT_EXPORTED_FROM_DLL
DISABLE_WARNING_BASE_NOT_EXPORTED_FROM_DLL
#endif

typedef std::string                String;
typedef std::vector<double>        VectorDouble; /// TODO : Create a class (fill, sum, mean...)
typedef std::vector<int>           VectorInt;
typedef std::vector<bool>          VectorBool;
typedef std::vector<float>         VectorFloat;
typedef std::vector<unsigned char> VectorUChar;
typedef std::vector<std::string>   VectorString;
typedef std::vector<VectorDouble>  VectorVectorDouble;
typedef std::vector<VectorInt>     VectorVectorInt;
typedef std::vector<VectorFloat>   VectorVectorFloat;

