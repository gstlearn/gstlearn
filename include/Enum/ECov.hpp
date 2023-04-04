/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
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
                 WENDLAND0,  24, "Wendland covariance (2,0)",\
                 WENDLAND1,  25, "Wendland covariance (3,1)",\
                 WENDLAND2,  26, "Wendland covariance (4,2)",\
                 MARKOV,     27, "Markovian covariances"

ENUM_DECLARE(ENUM_COV)
