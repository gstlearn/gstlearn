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
#include "Basic/Convolution.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Core/fftn.hpp"
#include "Basic/FFT.hpp"
#include "Db/DbGrid.hpp"
#include "Matrix/MatrixRectangular.hpp"

Convolution::Convolution(DbGrid* dbgrid)
  : _dbgrid(dbgrid)
{
}

Convolution::Convolution(const Convolution& m)
  : _dbgrid(m._dbgrid)
{
}

Convolution& Convolution::operator=(const Convolution &m)
{
  if (this != &m)
  {
    _dbgrid = m._dbgrid;
  }
  return *this;
}

Convolution::~Convolution()
{

}

bool Convolution::_isDbGridDefined() const
{
  if (_dbgrid == nullptr)
  {
    messerr("You must define 'dbgrid' beforehand");
    return false;
  }
  return true;
}

int Convolution::ConvolveSparse(int iatt,
                                const VectorVectorInt& ranks,
                                const MatrixRectangular& wgt,
                                const VectorDouble& means)
{
  if (! _isDbGridDefined()) return 1;
  int ndim    = _dbgrid->getNDim();
  int nvar    = _dbgrid->getNLoc(ELoc::Z);
  int nbneigh = (int)ranks.size();

  // Preliminary checks
  if (ndim != (int)ranks[0].size())
  {
    messerr("The second dimension of 'ranks' (%d)", (int)ranks[0].size());
    messerr("must be equal to the space dimension (%d)", ndim);
    return 1;
    }
  if (wgt.getNRows() != nbneigh * nvar)
  {
    messerr("The number of rows in the weight matrix (%d)", wgt.getNRows());
    messerr("must be equal to the number of neighbors (%d)", nbneigh);
    messerr("times the number of variables (%d)", nvar);
    return 1;
  }
  if (wgt.getNCols() != nvar)
  {
    messerr("The number of columns in the weight matrix (%d)", wgt.getNCols());
    messerr(" must be equal to the number of variables (%d)", nvar);
    return 1;
  }

  VectorInt indTarget(ndim);
  VectorInt current(ndim);
  VectorDouble data(nvar);
  VectorDouble result(nvar);

  // Loop on the target samples
  for (int iech = 0, nech = _dbgrid->getNSample(); iech < nech; iech++)
  {
    if (!_dbgrid->isActive(iech)) continue;
    _dbgrid->rankToIndice(iech, indTarget);

    // Loop on the neighborhing samples
    result.fill(0.);
    bool correct = true;
    for (int ineigh = 0; ineigh < nbneigh && correct; ineigh++)
    {
      VH::addInPlace(indTarget, ranks[ineigh], current);
      correct = _dbgrid->getGrid().isInside(current);

      // Target is not estimated when one neighborhood sample is out of grid
      if (!correct) continue;

      // Load the values at location 'jech' for all variables and skip
      bool valid = true;
      int jech   = _dbgrid->indiceToRank(current);
      for (int ivar = 0; ivar < nvar && valid; ivar++)
      {
        double value = _dbgrid->getZVariable(jech, ivar);
        if (FFFF(value))
          valid = false;
        else
        {
          if (!means.empty()) value -= means[ivar];
          data[ivar] = value;
        }
      }

      // Target is not estimated if one data value is undefined
      if (!valid)
      {
        correct = false;
        continue;
      }

      // Loop on data variable
      for (int jvar = 0; jvar < nvar; jvar++)
      {
        int irow = jvar * nbneigh + ineigh;
          
        // Loop on target variable
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          result[ivar] += data[jvar] * wgt.getValue(irow, ivar);
        }
      }
    }

    // Store the results
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      double value = TEST;
      if (correct)
      {
        value = result[ivar];
        if (!means.empty()) value += means[ivar];
      }
      _dbgrid->setArray(iech, iatt + ivar, value);
    }
  }
  return 0;
}

int Convolution::ConvolveFFT(int iatt,
                             int nvar,
                             const DbGrid* marpat,
                             const VectorDouble& means)
{
  int ndim       = _dbgrid->getNDim();
  int nv2        = nvar * nvar;
  VectorInt dims = _dbgrid->getNXs();
  int sImage     = VH::product(dims);
  double rImage  = (double) sImage;

  // Find the center of kernel
  VectorInt cKernel = marpat->getCenterIndices(false);

  // Find the center of the image
  VectorInt cImage = _dbgrid->getCenterIndices(false);

  // Find the shift between the two centers
  VectorInt shift = VH::subtract(cKernel, cImage);

  // For each kernel, allocate arrays (real and imaginary parts)
  // at the dimension of the final image. 
  VectorInt indices(ndim);
  std::vector<VectorDouble> kernelRe;
  std::vector<VectorDouble> kernelIm;
  kernelRe.resize(nv2);
  kernelIm.resize(nv2);
  for (int iv2 = 0; iv2 < nv2; iv2++)
  {
    kernelRe[iv2].fill(0., sImage);
    kernelIm[iv2].fill(0., sImage);
  }

  // For each pair or variable (ivar, jvar), plunge the corresponding
  // set of weights (given by 'marpat') in the "middle" of each 'kernel'
  int sKernel = VH::product(marpat->getNXs());
  for (int i = 0; i < sKernel; i++)
  {
    marpat->rankToIndice(i, indices);
    VH::addInPlace(indices, shift);
    int j = _dbgrid->indiceToRank(indices);

    for (int iv2 = 0; iv2 < nv2; iv2++)
      kernelRe[iv2][j] = marpat->getLocVariable(ELoc::Z, i, iv2);
  }

  // Perform the FFT forward transform for each kernel
  for (int iv2 = 0; iv2 < nv2; iv2++)
    if (fftn(ndim, dims.data(), kernelRe[iv2].data(), kernelIm[iv2].data(), 1)) return 1;
  
  // Retreive the vector of images (real and imaginary parts)
  std::vector<VectorDouble> imageRe;
  std::vector<VectorDouble> imageIm;
  imageRe.resize(nvar);
  imageIm.resize(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    imageRe[ivar] = _dbgrid->getColumnByLocator(ELoc::Z, ivar);
    if (!means.empty()) VH::addConstant(imageRe[ivar], -means[ivar]);
    imageIm[ivar].resize(sImage, 0.);

    // Perform the FFT forward transform of each image
    if (fftn(ndim, dims.data(), imageRe[ivar].data(), imageIm[ivar].data(), 1)) return 1;
  }

  VectorDouble result(sImage);
  VectorDouble localRe(sImage);
  VectorDouble localIm(sImage);

  // Loop on the output variables
  int iv2 = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    result.fill(0.);
    for (int jvar = 0; jvar < nvar; jvar++, iv2++)
    {

      // Perform the element-wise product of complex vectors
      VH::multiplyComplexInPlace(kernelRe[iv2], kernelIm[iv2], imageRe[jvar],
                                 imageIm[jvar], localRe, localIm);

      // Compute the inverse FFT
      if (fftn(ndim, dims.data(), localRe.data(), localIm.data(), -1, rImage)) return 1;

      // Cumulate the real parts
      VH::addInPlace(result, localRe);
    }

    // Perform the ultimate swap
    fftshift(dims, result);

    // Store the results in the Db
    if (!means.empty()) VH::addConstant(result, means[ivar]);
    _dbgrid->setArrayByUID(result, iatt + ivar);
  }
  return 0;
}
