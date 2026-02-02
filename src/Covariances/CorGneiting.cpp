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
#include "Basic/LawStable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CorGaussianMixture.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "geoslib_define.h"
#include <cmath>
#include <memory>

namespace gstlrn
{

CorGneiting::CorGneiting(const CovContext& ctxt, const ECov& type)
  : ACov(ctxt)
  , _corS()
  , _corT()
  , _cor(*_corS)
  , _separability(1.0) {
    DECLARE_UNUSED(type)
  };

CorGneiting::CorGneiting(const CorGaussianMixture* covS, const CorAniso* covT, double separability)
  : ACov()
  , _corS(std::shared_ptr<const CorGaussianMixture>(std::dynamic_pointer_cast<const CorGaussianMixture>((covS->cloneShared()))))
  , _corT(std::shared_ptr<const CorAniso>(std::dynamic_pointer_cast<const CorAniso>((covT->cloneShared()))))
  , _cor(*covS)
  , _separability(separability)
{
  if (_separability < 0.0 || _separability > 1.0)
  {
    _separability = 0.0;
    messerr("CorGneiting: Separability must be in [0,1]");
    messerr("It has been set to 0");
  }
  _corS->setOptimEnabled(false);
  _corT->setOptimEnabled(false);
  _cor.setOptimEnabled(false);

  // Define the Context
  auto space = SpaceComposite::create();
  space->addSpaceComponent(covS->getSpace());
  space->addSpaceComponent(covT->getSpace());
  _ctxt.setSpace(space);

  Id nvar = covS->getNVar();
  CovContext ctxt(nvar, space);
  setContext(ctxt);
}

CorGneiting::CorGneiting(const CorGneiting& r)
  : ACov(r)
  , _corS(r._corS)
  , _corT(r._corT)
  , _cor(r._cor)
  , _separability(r._separability)
{
}

CorGneiting& CorGneiting::operator=(const CorGneiting& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    _corS         = r._corS;
    _corT         = r._corT;
    _cor          = r._cor;
    _separability = r._separability;
  }
  return *this;
}

CorGneiting::~CorGneiting()
{
}

CorGneiting* CorGneiting::create(
  const CovContext& ctxt, // spatial context
  const ECov& type,
  double alpha,
  double beta,
  double timeRange,
  const VectorDouble& params, // Nus
  const VectorDouble& kappas, // Kappas
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

  Id ndim = ctxt.getNDim(); // dimension of the space
  Id nvar = ctxt.getNVar();
  if (nvar != params.length())
  {
    messerr("Inconsistent number of shape parameters: nvar = %d vs. params.length = %d", nvar, params.length());
    return nullptr;
  }
  if (nvar != kappas.length())
  {
    messerr("Inconsistent number of rate parameters: nvar = %d vs. coeffScales.length = %d", nvar, kappas.length());
    return nullptr;
  }

  if ((alpha <= 0.0) || (alpha > 2.0))
  {
    messerr("Gneiting-Gaussian model requires alpha in (0,2]:  alpha = %f", alpha);
    return nullptr;
  }
  if ((beta <= 0.0) || (beta > 1.0))
  {
    messerr("Gneiting-Gaussian model requires beta in (0,2]: beta = %f", beta);
    return nullptr;
  }
  // parameters of the spatial and time anisotropy
  if (ctxt.getNDim() < 1)
  {
    messerr("Space-time model is defined for ndim >= 1");
    return nullptr;
  }

  Id nr = ranges.length();
  if (nr > 0)
  {
    if (nr != ndim)
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
  auto* corS = new CorGaussianMixture(
    CovContext(nvar, ndim),
    ECov::fromValue(id_spatial),
    params,
    kappas,
    ranges,
    angles,
    flagRange);
  // creation of the temporal covariance
  auto* corT = new CorAniso(CovContext(1, 1), ECov::CAUCHY_GEN);
  corT->setParam(alpha, 0);           // alpha in (0,2]
  corT->setParam(beta * ndim / 2, 1); // beta*d/2 with beta in (0,1]
  if (flagRange)
  {
    corT->setRange(0, timeRange);
  }
  else
  {
    corT->setScaleDim(0, timeRange);
  }
  return new CorGneiting(corS, corT, separability);
}

SpectrumRN CorGneiting::simulateSpectrumRN(Id ns, const ACov* cov0) const
{
  // simulation of the spatial spectrum and variables weights
  SpectrumRN sp      = _corS->simulateSpectrumRN(ns, cov0);
  Id ndim            = _corS->getNDim(); // spatial dimension
  MatrixDense omegaS = sp.getOmega();
  MatrixDense gamma  = sp.getGamma();

  // simulation of the space-time spectrum

  MatrixDense omega(ns, ndim + 1);
  for (Id idim = 0; idim < ndim; idim++)
  {
    omega.setColumn(idim, omegaS.getColumn(idim));
  }

  // simulation of the time frequencies
  double alpha       = _corT->getCorFunc()->getParam(0);
  double beta        = _corT->getCorFunc()->getParam(1) * 2 / ndim; // the second parameter of the Cauchy is beta*ndim/2
  double timeRange   = _corT->getScale(0);
  MatrixDense omega0 = sp.getOmega0(); // spatial frequency without anisotropy
  VectorDouble xi0   = sp.getXi();     // random scales of the Gaussian mixture

  for (Id ib = 0; ib < ns; ib++)
  {
    double no2 = 0.0;
    double xi  = xi0[ib];
    for (Id idim = 0; idim < ndim; idim++)
    {
      double val = omega0.getValue(ib, idim);
      no2 += (val * val);
    }
    double lb    = std::pow(no2 / (4.0 * xi), 1 / beta); // = lambda ^ (1/beta)
    double sim_S = LawStable::law_stable_unilateral_exptilt(beta, lb);
    double sim_T = LawStable::law_stable_bilateral(alpha);
    double val   = std::pow(sim_S * lb, 1 / alpha) * sim_T / timeRange;
    omega.setValue(ib, ndim, val);
  }
  return SpectrumRN(gamma, omega);
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
  Id ndim      = _corS->getNDim();
  double ct    = _corT->evalCov(p1_T, p2_T, 0, 0, mode);
  double scale = pow(ct, _separability / ndim);
  _cor.setScaleGneiting(scale);
  /*
  for (Id idim = 0; idim < ndim; idim++)
  {
    double val = _corS->getScaleDim(idim);
    _cor.setScaleCor(idim, val / scale);
  }
  */
  double cs = _cor.evalCov(p1_S, p2_S, ivar, jvar, mode);
  return cs * ct;
}
} // namespace gstlrn
