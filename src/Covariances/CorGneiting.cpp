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
                                 const VectorDouble& params,
                                 const VectorDouble& ranges,
                                 const VectorDouble& angles,
                                 double separability,
                                 bool flagRange)
{
  DECLARE_UNUSED(type)
  if (ctxt.getNVar() != 1)
  {
    messerr("This function is dedicated to the Monovariate case");
    return nullptr;
  }

  if (ctxt.getNDim() <= 1)
  {
    messerr("Space-time model is defined for ndim >= 2");
    return nullptr;
  }
  Id ndim = ctxt.getNDim() - 1; // dimension of the space

  // parameters of the trace covariances
  if (params.length() != 3)
  {
    messerr("Gneiting model requires 3 parameters");
    return nullptr;
  }
  double nu    = params[0];
  double alpha = params[1];
  double beta  = params[2];

  // parameters of the spatial anisotropy and the time scale
  VectorDouble rangesS(ndim, 1.0);
  VectorDouble rangesT(1, 1.0);
  Id nr = ranges.length();
  if (nr > 0)
  {
    if (nr != ndim + 1)
    {
      messerr("Inconsistent number of ranges (%d)", nr);
      return nullptr;
    }
    for (Id idim = 0; idim < ndim; idim++)
    {
      rangesS[idim] = ranges[idim];
    }
    rangesT[0] = ranges[ndim];
  }

  VectorDouble anglesS(ndim, 1.0);
  Id nang = angles.length();
  if (nang > 0)
  {
    if (nang != ndim)
    {
      messerr("Inconsistent number of angles (%d)", nang);
      return nullptr;
    }
    for (Id idim = 0; idim < ndim; idim++)
    {
      anglesS[idim] = angles[idim];
    }
  }

  // separability coefficient (TODO: not used)
  if ((separability < 0.0) || (separability > 1.0))
  {
    messerr("Inconsistent separability coefficient = %f", separability);
    return nullptr;
  }

  // creation of the spatial covariance
  CovContext ctxt_S(1, ndim);
  CorAniso* corS = CorAniso::createAnisotropic(
    ctxt_S,
    ECov::MATERN,
    rangesS,
    nu,
    anglesS,
    flagRange);

  // creation of the temporal covariance
  CovContext ctxt_T(1, 1);
  CorAniso* corT = CorAniso::createAnisotropic(
    ctxt_S,
    ECov::CAUCHY_GEN,
    rangesT,
    beta,
    VectorDouble(),
    flagRange);

  return new CorGneiting(corS, corT, separability);
}

  MatrixDense CorGneiting::simulateSpectralOmega(Id nb) const
{
  Id ndim      = _corS->getNDim();
  double nu    = _corS->getParam();
  double alpha = 2.0;               // TODO: to be modified into _corT->getParam(0)
  double beta  = _corT->getParam(); // TODO: to be modified into _corT->getParam(1)
  MatrixDense omegaS(nb, ndim);     // spatial frequencies
  MatrixDense omegaT(nb, 1);        // temporal frequencies
  for (Id n = 0; n < nb; n++)
  {
    double Rval = 1.0;
    if (nu > 0)
    {
      Rval /= (4.0 * law_gamma(nu));
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
  auto p1_S = p1.spacePointOnSubspace(0);
  auto p2_S = p2.spacePointOnSubspace(0);
  auto p1_T = p1.spacePointOnSubspace(1);
  auto p2_T = p2.spacePointOnSubspace(1);
  double ct = _corT->evalCov(p1_T, p2_T, ivar, jvar, mode);

  double scale = pow(ct, _separability / _corSCopy.getNDim(0));
  for (Id i = 0; i < _corSCopy.getNDim(); i++)
    _corSCopy.setScaleDim(i, _corS->getScale(i) / scale);
  double cs = _corSCopy.evalCov(p1_S, p2_S, ivar, jvar, mode);

  return cs * ct;
}
} // namespace gstlrn
