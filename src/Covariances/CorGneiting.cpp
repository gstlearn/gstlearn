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
#include "Covariances/ACor.hpp"
#include "Covariances/CovProportional.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovContext.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "geoslib_define.h"
#include <vector>

CorGneiting::CorGneiting(const CorAniso* covS,const CorAniso* covTemp, double separability)
: ACor(CovContext())
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
  delete _space;
  SpaceComposite* space = new SpaceComposite();
  space->addSpaceComponent(covS->getSpace());
  space->addSpaceComponent(covTemp->getSpace()); 
  _space = space;
}


CorGneiting::CorGneiting(const CorGneiting& r):
ACor(r)
, _covS(r._covS)
, _covTemp(r._covTemp)
, _separability(r._separability)
, _covSCopy(*r._covS)
{
}

CorGneiting& CorGneiting::operator=(const CorGneiting &r)
{
  if (this != &r)
  {
    ACor::operator =(r);
    _ctxt = r._ctxt;
    _covS = r._covS;
    _covTemp = r._covTemp;
    _covSCopy = r._covSCopy;
    _separability = r._separability;
  }
  return *this;
}

CorGneiting::~CorGneiting()
{
}

void CorGneiting::_optimizationSetTarget(const SpacePoint &pt) const 
{
  DECLARE_UNUSED(pt)
  //_covS->optimizationSetTarget(pt.spacePointOnSubspace(0));
  //_covTemp->optimizationSetTarget(pt.spacePointOnSubspace(1));
}
  
void CorGneiting::optimizationSetTargetByIndex(int iech) const 
{
  DECLARE_UNUSED(iech)
  //_covS->optimizationSetTargetByIndex(iech);
  //_covTemp->optimizationSetTargetByIndex(iech);
}

void CorGneiting::_optimizationPreProcess(const std::vector<SpacePoint>& p) const 
{
  DECLARE_UNUSED(p)
 // _covS->_optimizationPreProcess(p);
 // _covTemp->_optimizationPreProcess(p);
}

void CorGneiting::optimizationPostProcess() const
{
  //_covS->optimizationPostProcess();
  //_covTemp->optimizationPostProcess();
}

double CorGneiting::eval(const SpacePoint& p1,
                    const SpacePoint& p2,
                    int ivar,
                    int jvar,
                    const CovCalcMode* mode) const
{
  auto p1_0 = p1.spacePointOnSubspace(0);
  auto p2_0 = p2.spacePointOnSubspace(0);
  auto p1_1 = p1.spacePointOnSubspace(1);
  auto p2_1 = p2.spacePointOnSubspace(1);
  double ct = _covTemp->eval(p1_1, p2_1, ivar, jvar, mode);

  double scale = pow(ct,_separability/_covSCopy.getNDim(0));
  for (int i = 0; i < (int) _covSCopy.getNDim(); i++)
  {
    _covSCopy.setScale(i, _covS->getScale(i) / scale);
  }
  double cs = _covSCopy.eval(p1_0, p2_0, ivar, jvar, mode);
  
  return cs * ct; 
}
