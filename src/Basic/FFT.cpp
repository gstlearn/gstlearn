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
#include "Basic/FFT.hpp"
#include "Arrays/Array.hpp"
#include "Basic/VectorNumT.hpp"
#include "Core/fftn.hpp"

#include <cmath>

namespace gstlrn
{
/****************************************************************************/
/*!
 **  Calculate the FFT in a space of dimension N
 **
 ** \return  Error return code
 **
 *****************************************************************************/
Id FFTn(Id ndim,
         const VectorInt& dims,
         VectorDouble& Re,
         VectorDouble& Im,
         Id iSign,
         double scaling)
{
  Id n = static_cast<Id>(Re.size());
  Im.resize(n, 0.);
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
                          Id N,
                          const std::function<std::complex<double>(VectorDouble, double)>& funcSpectrum)
{
  Id ndim = static_cast<Id>(hmax.size());
  VectorInt nxs(ndim);
  for (Id idim = 0; idim < ndim; idim++)
    nxs[idim] = N;
  Array array(nxs);

  Id ntotal = static_cast<Id>(pow(N, ndim));
  VectorDouble a(ndim);
  double coeff = 0;
  double prod  = 1.;

  for (Id idim = 0; idim < ndim; idim++)
  {
    coeff   = 1. / (2. * hmax[idim]);
    a[idim] = GV_PI * (N - 1) / (hmax[idim]);
    prod *= coeff;
  }

  VectorDouble Re(ntotal);
  VectorDouble Im(ntotal);
  VectorInt indices(ndim);
  VectorDouble temp(ndim);
  for (Id iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad, indices);

    Id s = 1;
    for (Id idim = 0; idim < ndim; idim++)
    {
      temp[idim] = a[idim] * (static_cast<double>(indices[idim]) / (N - 1) - 0.5);
      s *= (indices[idim] % 2) ? -1 : 1;
    }

    std::complex<double> fourier = funcSpectrum(temp, time);
    Re[iad]                      = s * prod * fourier.real();
    Im[iad]                      = s * prod * fourier.imag();
  }
  FFTn(ndim, nxs, Re, Im);

  // Retrieve information from the Re array and load them back in the array result.

  for (Id iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad, indices);
    Id s = 1;
    for (Id idim = 0; idim < ndim; idim++)
    {
      s *= (indices[idim] % 2) ? -1 : 1;
    }
    array.setValue(indices, Re[iad] * s);
  }
  return array;
}

Array evalCovFFTSpatial(const VectorDouble& hmax,
                        Id N,
                        const std::function<double(const VectorDouble&)>& funcSpectrum)
{
  Id ndim = static_cast<Id>(hmax.size());
  VectorInt nxs(ndim);
  for (Id idim = 0; idim < ndim; idim++)
    nxs[idim] = N;
  Array array(nxs);

  Id ntotal = static_cast<Id>(pow(N, ndim));
  VectorDouble a(ndim);
  double coeff = 0;
  double prod  = 1.;

  for (Id idim = 0; idim < ndim; idim++)
  {
    coeff   = 1. / (2. * hmax[idim]);
    a[idim] = GV_PI * (N - 1) / (hmax[idim]);
    prod *= coeff;
  }

  VectorDouble Re(ntotal);
  VectorDouble Im(ntotal, 0.);
  VectorInt indices(ndim);
  VectorDouble temp(ndim);

  for (Id iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad, indices);

    Id s = 1;
    for (Id idim = 0; idim < ndim; idim++)
    {
      temp[idim] = a[idim] * (static_cast<double>(indices[idim]) / (N - 1) - 0.5);
      s *= (indices[idim] % 2) ? -1 : 1;
    }
    Re[iad] = s * prod * funcSpectrum(temp);
  }
  FFTn(ndim, nxs, Re, Im);

  // Retrieve information from the Re array and load them back in the array result.

  for (Id iad = 0; iad < ntotal; iad++)
  {
    array.rankToIndice(iad, indices);
    Id s = 1;
    for (Id idim = 0; idim < ndim; idim++)
    {
      s *= (indices[idim] % 2) ? -1 : 1;
    }
    array.setValue(indices, Re[iad] * s);
  }

  return array;
}

// Fonction pour calculer les strides d'un tableau en mémoire contiguë (C row-major)
static VectorInt _computeStrides(Id ndim, const VectorInt& dims)
{
  VectorInt strides(ndim);
  Id stride = 1;
  for (Id idim = 0; idim < ndim; idim++)
  {
    strides[idim] = stride;
    stride *= dims[idim];
  }
  return strides;
}

// Fonction pour calculer les tailles des demi-champs
static VectorInt _computeHalf(Id ndim, const VectorInt& dims)
{
  VectorInt half(ndim);
  for (Id idim = 0; idim < ndim; idim++)
  {
    half[idim] = ceil(dims[idim] / 2.);
  }
  return half;
}

static Id _getIndex(Id ndim, const VectorInt& strides, const VectorInt& indices)
{
  Id index = 0;
  for (Id idim = 0; idim < ndim; idim++)
    index += indices[idim] * strides[idim];
  return index;
}

// Fonction pour effectuer un fftshift sur un tenseur nD stocké en 1D
void fftshift(const VectorInt& dims, VectorDouble& data)
{
  Id ndim       = static_cast<Id>(dims.size());
  Id total_size = static_cast<Id>(data.size());

  // Calcul des moitiés des dimensions
  VectorInt half = _computeHalf(ndim, dims);

  // Calcule des strides
  VectorInt strides = _computeStrides(ndim, dims);

  // Réorganisation des indices : swap des quadrants
  VectorDouble temp(data);
  VectorInt coords_old(ndim);
  VectorInt coords_new(ndim);
  for (Id index = 0; index < total_size; index++)
  {
    // Convertir index linéaire en coordonnées nD
    Id linear_index = index;
    for (Id idim = 0; idim < ndim; idim++)
    {
      coords_old[idim] = linear_index % dims[idim];
      linear_index /= dims[idim];
    }

    // Appliquer fftshift en inversant les moitiés
    for (Id dim = 0; dim < ndim; dim++)
      coords_new[dim] = (coords_old[dim] + half[dim]) % dims[dim];

    // Reconvertir en index linéaire
    Id new_index   = _getIndex(ndim, strides, coords_new);
    data[new_index] = temp[index];
  }
}
} // namespace gstlrn
