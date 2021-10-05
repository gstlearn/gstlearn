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

#include "Enum/AEnum.hpp"

#define ENUM_COV ECov, UNKNOWN,\
                 UNKNOWN,    -2, "Unknown covariance",\
                 FUNCTION,   -1, "External covariance function",\
                 NUGGET,      0, "Nugget effect",\
                 EXPONENTIAL, 1, "Exponential",\
                 SPHERICAL,   2, "Spherical",\
                 GAUSSIAN,    3, "Gaussian",\
                 CUBIC,       4, "Cubic",\
                 SINCARD,     5, "Sine Cardinal",\
                 BESSEL_J,    6, "Bessel J",\
                 BESSEL_K,    7, "Bessel K",\
                 GAMMA,       8, "Gamma",\
                 CAUCHY,      9, "Cauchy",\
                 STABLE,     10, "Stable",\
                 LINEAR,     11, "Linear",\
                 POWER,      12, "Power",\
                 ORDER1_GC,  13, "First Order Generalized covariance",\
                 SPLINE_GC,  14, "Spline Generalized covariance",\
                 ORDER3_GC,  15, "Third Order Generalized covariance",\
                 ORDER5_GC,  16, "Fifth Order Generalized covariance",\
                 COSINUS,    17, "Cosine",\
                 TRIANGLE,   18, "Triangle",\
                 COSEXP,     19, "Cosine Exponential",\
                 REG1D,      20, "1-D Regular",\
                 PENTA,      21, "Pentamodel",\
                 SPLINE2_GC, 22, "Order-2 Spline",\
                 STORKEY,    23, "Storkey covariance in 1-D",\
                 WENDLAND1,  24, "Wendland covariance (first type)",\
                 WENDLAND2,  25, "Wendland covariance (second type)",\
                 P8,         26, "Polynomial of degree 8"

ENUM_DECLARE(ENUM_COV)
