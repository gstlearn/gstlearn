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
#include "geoslib_old_f.h"

#include "Basic/AException.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/FFT.hpp"
#include "Arrays/Array.hpp"

#include <math.h>

/****************************************************************************/
/*!
 **  Calculate the FFT in a space of dimension N
 **
 ** \return  Error return code
 **
 *****************************************************************************/
int FFTn(int ndim,
         const VectorInt& dims,
         VectorDouble& Re,
         VectorDouble& Im,
         int iSign,
         double scaling)
{
  int n = (int) Re.size();
  Im.resize(n,0.);
  return fftn(ndim, dims.data(), Re.data(), Im.data(), iSign, scaling);
}

/**
 * perform the FFT transform for a First-Order Space Time evolution equation
 * @param hmax Maximum spatial distances (Dimension: spatial ndim)
 * @param time Time of the covariance slice
 * @param N    Discretization number (in each spatial dimension)
 * @param funcSpectrum External adequate spectrum evaluation function
 * @return Array of spatio-temporal covariance
 */
Array evalCovFFTTimeSlice(const VectorDouble& hmax, double time, int N,
                          std::function<std::complex<double>(VectorDouble, double)> funcSpectrum)
{
  int ndim = (int) hmax.size();
  VectorInt nxs(ndim);
  for (int idim = 0; idim < ndim; idim++)
    nxs[idim] = N ;
  Array array(nxs);

  int ntotal = (int) pow(N, ndim);
  VectorDouble a(ndim);
  double coeff = 0;
  double prod = 1.;

  for(int idim = 0; idim < ndim; idim++)
  {
    coeff = 1. / (2. * hmax[idim]);
    a[idim] = GV_PI * (N-1) / (hmax[idim]);
    prod *= coeff;
  }

  VectorDouble Re(ntotal);
  VectorDouble Im(ntotal);
  VectorInt indices(ndim);
  VectorDouble temp(ndim);
  for (int iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad,indices);

    int s = 1;
    for (int idim = 0; idim < ndim; idim++)
    {
      temp[idim] = a[idim] * ((double)indices[idim] / (N - 1) - 0.5);
      s *= (indices[idim] % 2) ? -1 : 1;
    }

    std::complex<double> fourier = funcSpectrum(temp,time);
    Re[iad] = s * prod * fourier.real();
    Im[iad] = s * prod * fourier.imag();
  }
  FFTn(ndim, nxs, Re, Im);

  // Retrieve information from the Re array and load them back in the array result.

  for (int iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad,indices);
    int s = 1;
    for (int idim = 0;  idim < ndim; idim++)
    {
      s *= (indices[idim] % 2) ? -1 : 1;
    }
    array.setValue(indices,Re[iad] * s);
  }
  return array;
}

Array evalCovFFTSpatial(const VectorDouble& hmax, int N,
                        std::function<double(const VectorDouble&)> funcSpectrum)
{
  int ndim = (int) hmax.size();
  VectorInt nxs(ndim);
  for (int idim = 0; idim < ndim; idim++)
    nxs[idim] = N;
  Array array(nxs);

  int ntotal = (int) pow(N, ndim);
  VectorDouble a(ndim);
  double coeff = 0;
  double prod = 1.;

  for (int idim = 0; idim < ndim; idim++)
  {
    coeff = 1. / (2. * hmax[idim]);
    a[idim] = GV_PI * (N - 1) / (hmax[idim]);
    prod *= coeff;
  }

  VectorDouble Re(ntotal);
  VectorDouble Im(ntotal, 0.);
  VectorInt indices(ndim);
  VectorDouble temp(ndim);

  for (int iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad, indices);

    int s = 1;
    for (int idim = 0; idim < ndim; idim++)
    {
      temp[idim] = a[idim] * ((double) indices[idim] / (N - 1) - 0.5);
      s *= (indices[idim] % 2) ? -1 : 1;
    }
    Re[iad] = s * prod * funcSpectrum(temp);
  }
  FFTn(ndim, nxs, Re, Im);

  // Retrieve information from the Re array and load them back in the array result.

  for (int iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad,indices);
    int s = 1;
    for (int idim = 0;  idim < ndim; idim++)
    {
      s *= (indices[idim] % 2) ? -1 : 1;
    }
    array.setValue(indices,Re[iad] * s);
  }

  return array;
}
