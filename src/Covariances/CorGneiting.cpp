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

#include "Covariances/CorGneiting.hpp"
#include "Basic/AStringable.hpp"
#include "Covariances/CovContext.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "Covariances/CovCalcMode.hpp"
#include <vector>

CorGneiting::CorGneiting(const CorAniso* covS, const CorAniso* covTemp, double separability)
  : ACov()
  , _covS(covS)
  , _covTemp(covTemp)
  , _separability(separability)
  , _covSCopy(*covS)
{
  if (separability < 0.0 || separability > 1.0)
  {
    _separability = 0;
    messerr("CorGneiting: Separability must be in [0,1]");
    messerr("It has been set to 0");
  }
  _covS->setOptimEnabled(false);
  _covSCopy.setOptimEnabled(false);
  _covTemp->setOptimEnabled(false);

  // Define the Context
  auto space = SpaceComposite::create();
  space->addSpaceComponent(covS->getSpace());
  space->addSpaceComponent(covTemp->getSpace());
  _space = space;

  int nvar = covS->getNVar();
  CovContext ctxt    = CovContext(nvar, space);
  setContext(ctxt);
}

CorGneiting::CorGneiting(const CorGneiting& r)
  : ACov(r)
  , _covS(r._covS)
  , _covTemp(r._covTemp)
  , _separability(r._separability)
  , _covSCopy(*r._covS)
{
}

CorGneiting& CorGneiting::operator=(const CorGneiting& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    _covS         = r._covS;
    _covTemp      = r._covTemp;
    _covSCopy     = r._covSCopy;
    _separability = r._separability;
  }
  return *this;
}

CorGneiting::~CorGneiting()
{
}

// void CorGneiting::_optimizationSetTarget(SpacePoint& pt) const
// {
//   _covS->optimizationSetTarget(pt);
//   _covSCopy.optimizationSetTarget(pt);
//   _covTemp->optimizationSetTarget(pt);
// }

// void CorGneiting::_optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const
// {
//   // DECLARE_UNUSED(mode)
//   // DECLARE_UNUSED(ps)
//   _covS->optimizationPreProcess(mode, ps);
//   _covSCopy.optimizationPreProcess(mode, ps);
//   _covTemp->optimizationPreProcess(mode, ps);
// }

// void CorGneiting::_optimizationPostProcess() const
// {
//   _covS->optimizationPostProcess();
//   _covSCopy.optimizationPostProcess();
//   _covTemp->optimizationPostProcess();
// }

double CorGneiting::_eval(const SpacePoint& p1,
                          const SpacePoint& p2,
                          int ivar,
                          int jvar,
                          const CovCalcMode* mode) const
{
  auto p1_S = p1.spacePointOnSubspace(0);
  auto p2_S = p2.spacePointOnSubspace(0);
  auto p1_T = p1.spacePointOnSubspace(1);
  auto p2_T = p2.spacePointOnSubspace(1);
  double ct = _covTemp->evalCov(p1_T, p2_T, ivar, jvar, mode);

  double scale = pow(ct, _separability / _covSCopy.getNDim(0));
  for (int i = 0; i < (int)_covSCopy.getNDim(); i++)
    _covSCopy.setScale(i, _covS->getScale(i) / scale);
  double cs = _covSCopy.evalCov(p1_S, p2_S, ivar, jvar, mode);

  return cs * ct;
}
