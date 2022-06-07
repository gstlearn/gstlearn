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
#include "Basic/Vector.hpp"
#include "Basic/Array.hpp"

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
