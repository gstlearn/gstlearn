/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Simulation/CalcSimuTurningBands.hpp"
#include "Covariances/CovHelper.hpp"
#include "Covariances/ACovFunc.hpp"
#include "Covariances/CovFactory.hpp"
#include "Basic/VectorT.hpp"
#include "Enum/ECov.hpp"

bool _isSelected(ACovFunc* cov, int ndim, int minorder, bool hasrange, bool flagSimtub)
{
  if (cov == nullptr) return false;
  if (ndim > (int) cov->getMaxNDim()) return false;
  if (minorder < cov->getMinOrder()) return false;
  if (hasrange && ! cov->hasRange()) return false;
  if (flagSimtub && ! isCovValidForTurningBands(cov->getType())) return false;
  return true;
}

VectorString CovHelper::getAllCovariances(int ndim, int minorder, bool hasrange, bool flagSimtub)
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
    if (_isSelected(cov, ndim, minorder, hasrange, flagSimtub))
      vs.push_back(it.getKey());

    delete cov;
    it.toNext();
  }
  return vs;
}
