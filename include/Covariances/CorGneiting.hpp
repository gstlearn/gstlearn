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
#pragma once

#include "Basic/ICloneable.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CorGaussianMixture.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class ACov;
class CorAniso;
/**
 * \brief
 * This class describes the univariate Gneiting correlation function.
 *
 */

class GSTLEARN_EXPORT CorGneiting: public ACov
{
public:
  CorGneiting(const CovContext& ctxt, const ECov& type);
  CorGneiting(const CorGaussianMixture* covS, const CorAniso* covT, double separability = 1.0);
  CorGneiting(const CorGneiting& r);
  CorGneiting& operator=(const CorGneiting& r);
  virtual ~CorGneiting();

  static CorGneiting* create(
    const CovContext& ctxt,
    const ECov& type,
    double alpha,
    double beta,
    double timeRange,           // time scale or range according to flagRange value
    const VectorDouble& params, // shape parameter for each variable
    const VectorDouble& kappas, // scale parameter for each variable
    const VectorDouble& ranges, // spatial scales or ranges according to flagRange value
    const VectorDouble& angles = VectorDouble(),
    double separability        = 1.0, // in [0,1] with a time-space separable model for 0, and "pure" Gneiting for 1
    bool flagRange             = false);

  IMPLEMENT_CLONING(CorGneiting)

  /// ACov Interface
  bool isConsistent(const ASpace* space) const override
  {
    return (space->getType() == ESpaceType::RN);
  }

  Id getNFac() const override;
  Id getNVar() const override { return _corS->getNVar(); }
  double getC0(Id ivar, Id jvar) const { return _corS->getC0(ivar, jvar); }
  double getNu(Id ivar, Id jvar) const { return _corS->getNu(ivar, jvar); }
  double getKappa(Id ivar, Id jvar) const { return _corS->getKappa(ivar, jvar); }
  double separability() const { return _separability; }

  bool isValidForSimulation(const ESimuType& simuType) const override
  {
    return (simuType == ESimuType::SPECTRAL);
  }
  SpectrumRN simulateSpectrumRN(Id ns, const ACov* cov0 = nullptr) const override;
  SpectrumOnRN* simulateOnRN(Id ns) const override;

  SpectrumOnRN* simulateSpace(Id ns) const {return _corS->simulateOnRN(ns);};
  SpectrumOnRN* simulateTime(Id ns) const {return _corT->simulateOnRN(ns);};

protected:
  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;

private:
  std::shared_ptr<const CorGaussianMixture> _corS; // the space covariance function (ndim, nvar >=1)
  std::shared_ptr<const CorAniso> _corT;           // the time covariance function (ndim = 1, nvar = 1)
  mutable CorGaussianMixture _cor;                 // an auxiliary space covariance function
  double _separability;
};

} // namespace gstlrn