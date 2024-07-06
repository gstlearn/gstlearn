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
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Covariances/CovHelper.hpp"
#include "Covariances/ACovFunc.hpp"
#include "Covariances/CovFactory.hpp"
#include "Basic/VectorT.hpp"
#include "Enum/ECov.hpp"

bool _isSelected(ACovFunc *cov,
                 int ndim,
                 int minorder,
                 bool hasrange,
                 bool flagSimtub,
                 bool flagSimuSpectral)
{
  if (cov == nullptr) return false;
  if (ndim > (int) cov->getMaxNDim()) return false;
  if (minorder < cov->getMinOrder()) return false;
  if (hasrange && ! cov->hasRange()) return false;
  if (flagSimtub && ! cov->isValidForTurningBand()) return false;
  if (flagSimuSpectral && ! cov->isValidForSpectral()) return false;
  return true;
}

/**
 * Returns the list of covariance names that are valid according to various options:
 *
 * @param ndim     Space dimension
 * @param minorder Minimum degree of the IRF (-1 for stationary, 0 for Intrinsic, ...)
 * @param hasrange Check if the Covariance has a Range defined
 * @param flagSimtub Check that the Covariance can be simulated using Turning Band Method
 * @param flagSimuSpectral Check if the Covariance can be simulated using the Spectral Method
 */
VectorString CovHelper::getAllCovariances(int  ndim,
                                          int  minorder,
                                          bool hasrange,
                                          bool flagSimtub,
                                          bool flagSimuSpectral)
{
  VectorString vs;

  // Create the Context (using Space Dimension)
  const CovContext ctxt = CovContext(1,ndim);

  // Loop on the basic structures
  auto it = ECov::getIterator();
  while (it.hasNext())
  {
    const ECov& covType = ECov::fromKey(it.getKey());
    ACovFunc* cov = CovFactory::createCovFunc(covType, ctxt);

    // Check if the covariance is valid
    if (_isSelected(cov, ndim, minorder, hasrange, flagSimtub, flagSimuSpectral))
      vs.push_back(it.getKey());

    delete cov;
    it.toNext();
  }
  return vs;
}
