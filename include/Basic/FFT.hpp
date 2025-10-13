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

#include "Arrays/Array.hpp"
#include "gstlearn_export.hpp"

#include <cmath>
#include <complex>
#include <functional>

namespace gstlrn
{
GSTLEARN_EXPORT int FFTn(int ndim,
                         const VectorInt& dims,
                         VectorDouble& Re,
                         VectorDouble& Im,
                         int iSign      = 1,
                         double scaling = 1.);
GSTLEARN_EXPORT Array evalCovFFTTimeSlice(const VectorDouble& hmax, double time, int N, const std::function<std::complex<double>(VectorDouble, double)>& funcSpectrum);
GSTLEARN_EXPORT Array evalCovFFTSpatial(const VectorDouble& hmax, int N, const std::function<double(const VectorDouble&)>& funcSpectrum);
GSTLEARN_EXPORT void fftshift(const VectorInt& dims, VectorDouble& data);
} // namespace gstlrn