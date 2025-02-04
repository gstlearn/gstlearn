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
#include "Basic/VectorNumT.hpp"
#include "Basic/FFT.hpp"
#include "Arrays/Array.hpp"
#include "Core/fftn.hpp"

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
Array evalCovFFTTimeSlice(const VectorDouble& hmax,
                          double time,
                          int N,
                          const std::function<std::complex<double>(VectorDouble, double)>& funcSpectrum)
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

Array evalCovFFTSpatial(const VectorDouble &hmax,
                        int N,
                        const std::function<double(const VectorDouble&)>& funcSpectrum)
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

// Fonction pour calculer les strides d'un tableau en mémoire contiguë (C row-major)
static VectorInt _computeStrides(int ndim, const VectorInt& dims)
{
  VectorInt strides(ndim);
  int stride = 1;
  for (int idim = 0; idim < ndim; idim++)
  {
    strides[idim] = stride;
    stride *= dims[idim];
  }
  return strides;
}

// Fonction pour calculer les tailles des demi-champs
static VectorInt _computeHalf(int ndim, const VectorInt& dims)
{
  VectorInt half(ndim);
  for (int idim = 0; idim < ndim; idim++)
  {
    half[idim] = ceil(dims[idim] / 2.);
  }
  return half;
}

static int _getIndex(int ndim, const VectorInt& strides, const VectorInt& indices)
{
  int index = 0;
  for (int idim = 0; idim < ndim; idim++)
    index += indices[idim] * strides[idim];
  return index;
}

// Fonction pour effectuer un fftshift sur un tenseur nD stocké en 1D
void fftshift(const VectorInt& dims, VectorDouble& data)
{
  int ndim = (int) dims.size();
  int total_size = (int)data.size();

  // Calcul des moitiés des dimensions
  VectorInt half = _computeHalf(ndim, dims);

  // Calcule des strides
  VectorInt strides = _computeStrides(ndim, dims);

  // Réorganisation des indices : swap des quadrants
  VectorDouble temp(data);
  VectorInt coords_old(ndim);
  VectorInt coords_new(ndim);
  for (int index = 0; index < total_size; index++)
  {
    // Convertir index linéaire en coordonnées nD
    int linear_index = index;
    for (int idim = 0; idim < ndim; idim++)
    {
      coords_old[idim] = linear_index % dims[idim];
      linear_index /= dims[idim];
    }

    // Appliquer fftshift en inversant les moitiés
    for (int dim = 0; dim < ndim; dim++)
      coords_new[dim] = (coords_old[dim] + half[dim]) % dims[dim];

    // Reconvertir en index linéaire
    int new_index = _getIndex(ndim, strides, coords_new);
    data[new_index] = temp[index];
  }
}

