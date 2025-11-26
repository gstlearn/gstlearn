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
#include "Basic/Law.hpp"
#include "Basic/LawStable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include <cmath>
#include <memory>

namespace gstlrn
{

CorGneiting::CorGneiting(const ECov& type, const CovContext& ctxt)
  : ACov(ctxt)
  , _corS()
  , _corT()
  , _separability(0.0)
  , _corSCopy(*_corS) {
    DECLARE_UNUSED(type)
  };

CorGneiting::CorGneiting(const CorAniso* covS, const CorAniso* covTemp, double separability)
  : ACov()
  , _corS(std::shared_ptr<const CorAniso>(std::dynamic_pointer_cast<const CorAniso>((covS->cloneShared()))))
  , _corT(std::shared_ptr<const CorAniso>(std::dynamic_pointer_cast<const CorAniso>((covTemp->cloneShared()))))
  , _separability(separability)
  , _corSCopy(*covS)
{
  if (separability < 0.0 || separability > 1.0)
  {
    _separability = 0;
    messerr("CorGneiting: Separability must be in [0,1]");
    messerr("It has been set to 0");
  }
  _corS->setOptimEnabled(false);
  _corSCopy.setOptimEnabled(false);
  _corT->setOptimEnabled(false);

  // Define the Context
  auto space = SpaceComposite::create();
  space->addSpaceComponent(covS->getSpace());
  space->addSpaceComponent(covTemp->getSpace());
  _ctxt.setSpace(space);

  Id nvar = covS->getNVar();
  CovContext ctxt(nvar, space);
  setContext(ctxt);
}

CorGneiting::CorGneiting(const CorGneiting& r)
  : ACov(r)
  , _corS(r._corS)
  , _corT(r._corT)
  , _separability(r._separability)
  , _corSCopy(*r._corS)
{
}

CorGneiting& CorGneiting::operator=(const CorGneiting& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    _corS         = r._corS;
    _corT         = r._corT;
    _corSCopy     = r._corSCopy;
    _separability = r._separability;
  }
  return *this;
}

CorGneiting::~CorGneiting()
{
}

CorGneiting* CorGneiting::create(const ECov& type,
                                 const CovContext& ctxt,
                                 const VectorDouble& params, // = (alpha, beta*d/2, nu)
                                 const VectorDouble& ranges,
                                 const VectorDouble& angles,
                                 double separability,
                                 bool flagRange)
{
  Id id_spatial = 0;
  if ((type == ECov::GNEITING_G) || (type == ECov::GNEITING_C) || (type == ECov::GNEITING_M))
  {
    if (type == ECov::GNEITING_G)
    {
      id_spatial = ECov::GAUSSIAN.getValue();
    }
    if (type == ECov::GNEITING_M)
    {
      id_spatial = ECov::MATERN.getValue();
    }
    if (type == ECov::GNEITING_C)
    {
      id_spatial = ECov::CAUCHY.getValue();
    }
  }
  else
  {
    messerr("This function implements the Gneiting'model");
    return nullptr;
  }

  if (ctxt.getNVar() != 1)
  {
    messerr("This function is dedicated to the Monovariate case");
    return nullptr;
  }
  if (type == ECov::GNEITING_G)
  {
    if (params.length() != 2)
    {
      messerr("Gneiting-Gaussian model requires 2 parameters");
      return nullptr;
    }
  }
  else
  {
    if (params.length() != 3)
    {
      messerr("Gneiting model requires 3 parameters");
      return nullptr;
    }
  }
  // parameters of the spatial and time anisotropy
  if (ctxt.getNDim() < 2)
  {
    messerr("Space-time model is defined for ndim >= 2");
    return nullptr;
  }
  Id ndim = ctxt.getNDim() - 1; // dimension of the space
  Id nr   = ranges.length();
  if (nr > 0)
  {
    if (nr != ndim + 1)
    {
      messerr("Inconsistent number of ranges (%d vs. ndim = %d)", nr, ndim);
      return nullptr;
    }
  }
  Id nang = angles.length();
  if (nang > 0)
  {
    if (nang != ndim)
    {
      messerr("Inconsistent number of angles (%d)", nang);
      return nullptr;
    }
  }
  // separability coefficient
  if ((separability < 0.0) || (separability > 1.0))
  {
    messerr("Inconsistent separability coefficient = %f", separability);
    return nullptr;
  }

  // creation of the spatial covariance
  auto* corS = new CorAniso(ECov::fromValue(id_spatial), CovContext(1, ndim));
  if (type != ECov::GNEITING_G)
  {
    corS->setParam(params[2], 0); // nu > 0
  }
  for (Id idim = 0; idim < ndim; idim++)
  {
    if (flagRange)
    {
      corS->setRange(idim, ranges[idim]);
    }
    else
    {
      corS->setScaleDim(idim, ranges[idim]);
    }
  }
  if (nang > 0)
  {
    corS->setAnisoAngles(angles);
  }

  // creation of the temporal covariance
  auto* corT = new CorAniso(ECov::CAUCHY_GEN, CovContext(1, 1));
  corT->setParam(params[0], 0); // alpha in (0,2]
  corT->setParam(params[1], 1); // beta*d/2 with beta in (0,1]
  if (flagRange)
  {
    corT->setRange(0, ranges[ndim]);
  }
  else
  {
    corT->setScaleDim(0, ranges[ndim]);
  }
  return new CorGneiting(corS, corT, separability);
}

bool CorGneiting::isValidForSpectral() const
{
  return true;
}

MatrixDense CorGneiting::simulateSpectralOmega(Id nb) const
{
  Id ndim      = _corS->getNDim();
  double alpha = _corT->getCorFunc()->getParam(0);
  double beta  = _corT->getCorFunc()->getParam(1) * 2 / ndim;
  ECov type    = _corS->getCorFunc()->getType();
  double nu    = 0.0;
  if ((type == ECov::MATERN) || (type == ECov::CAUCHY))
  {
    nu = _corS->getCorFunc()->getParam(0);
  }
  MatrixDense omegaS(nb, ndim); // spatial frequencies
  MatrixDense omegaT(nb, 1);    // temporal frequencies
  for (Id n = 0; n < nb; n++)
  {
    double Rval = 1.0;
    if (type == ECov::MATERN)
    {
      Rval /= (4.0 * law_gamma(nu));
    }
    else if (type == ECov::CAUCHY)
    {
      Rval *= law_gamma(nu);
    }
    double s2r = std::sqrt(2.0 * Rval);
    double no2 = 0.0;
    for (Id idim = 0; idim < ndim; idim++)
    {
      double u = s2r * law_gaussian();
      omegaS.setValue(n, idim, u);
      no2 += u * u;
    }
    double lb    = std::pow(no2 / (4.0 * Rval), 1 / beta); // = lambda ^ (1/beta)
    double sim_S = LawStable::law_stable_unilateral_exptilt(beta, lb);
    double sim_T = LawStable::law_stable_bilateral(alpha);
    double val   = std::pow(sim_S * lb, 1 / alpha) * sim_T;
    omegaT.setValue(n, 0, val);
  }
  // applying the space anisotropy
  const auto& tensorS = _corS->getAniso().getTensorInverse();
  omegaS.prodMat(&tensorS);
  // applying the time anisotropy
  const auto& tensorT = _corT->getAniso().getTensorInverse();
  omegaT.prodMat(&tensorT);
  MatrixDense omegaST(nb, ndim + 1);
  for (Id icol = 0; icol < ndim; icol++)
  {
    omegaST.setColumn(icol, omegaS.getColumn(icol));
  }
  omegaST.setColumn(ndim, omegaT.getColumn(0));
  return omegaST;
}

// void CorGneiting::_optimizationSetTarget(SpacePoint& pt) const
// {
//   _corS->optimizationSetTarget(pt);
//   _corSCopy.optimizationSetTarget(pt);
//   _corT->optimizationSetTarget(pt);
// }

// void CorGneiting::_optimizationPreProcess(Id mode, const std::vector<SpacePoint>& ps) const
// {
//   // DECLARE_UNUSED(mode)
//   // DECLARE_UNUSED(ps)
//   _corS->optimizationPreProcess(mode, ps);
//   _corSCopy.optimizationPreProcess(mode, ps);
//   _corT->optimizationPreProcess(mode, ps);
// }

// void CorGneiting::_optimizationPostProcess() const
// {
//   _corS->optimizationPostProcess();
//   _corSCopy.optimizationPostProcess();
//   _corT->optimizationPostProcess();
// }

double CorGneiting::_eval(const SpacePoint& p1,
                          const SpacePoint& p2,
                          Id ivar,
                          Id jvar,
                          const CovCalcMode* mode) const
{
  auto p1_S    = p1.spacePointOnSubspace(0);
  auto p2_S    = p2.spacePointOnSubspace(0);
  auto p1_T    = p1.spacePointOnSubspace(1);
  auto p2_T    = p2.spacePointOnSubspace(1);
  double ct    = _corT->evalCov(p1_T, p2_T, ivar, jvar, mode);
  double scale = pow(ct, _separability / _corSCopy.getNDim(0));
  for (Id i = 0; i < _corSCopy.getNDim(); i++)
    _corSCopy.setScaleDim(i, _corS->getScale(i) / scale);
  double cs = _corSCopy.evalCov(p1_S, p2_S, ivar, jvar, mode);
  return cs * ct;
}
} // namespace gstlrn
