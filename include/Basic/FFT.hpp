/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Arrays/Array.hpp"

#include <math.h>
#include <complex>
#include <functional>

GSTLEARN_EXPORT int FFTn(int ndim,
                         const VectorInt& dims,
                         VectorDouble& Re,
                         VectorDouble& Im,
                         int iSign = 1,
                         double scaling = 1.);
GSTLEARN_EXPORT Array evalCovFFTTimeSlice(const VectorDouble& hmax, double time, int N,
                                          std::function<std::complex<double>(VectorDouble, double)> funcSpectrum);
GSTLEARN_EXPORT Array evalCovFFTSpatial(const VectorDouble& hmax, int N,
                                        std::function<double(const VectorDouble&)> funcSpectrum);
