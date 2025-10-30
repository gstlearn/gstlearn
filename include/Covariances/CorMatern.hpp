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
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include <memory>
#include <vector>

namespace gstlrn
{
class CorAniso;
/**
 * \brief
 * This class describes the multivariate Matern correlation function.
 *
 */
class GSTLEARN_EXPORT CorMatern: public ACov
{
public:
  CorMatern(const VectorDouble& ranges      = VectorDouble(),
            const VectorDouble& angles      = VectorDouble(),
            const VectorDouble& coeffScales = VectorDouble(),
            const VectorDouble& params      = VectorDouble(),
            const MatrixSymmetric& sigma    = MatrixSymmetric(),
            bool flagRange                  = true);
  CorMatern(const CorMatern& r);
  CorMatern& operator=(const CorMatern& r);
  virtual ~CorMatern();
  IMPLEMENT_CLONING(CorMatern)

  bool isConsistent(const ASpace* space) const override
  {
    DECLARE_UNUSED(space)
    return true;
  }
  bool isValidForSpectral() const override
  {
    return true;
  }

  const CorAniso* getCorRef() const { return _corRef.get(); }
  /// ACov Interface

  Id getNVar() const override { return _nVar; }
  double getC0(Id ivar, Id jvar) const { return _C0.getValue(ivar, jvar); }
  double computeScale(Id ivar, Id jvar) const;
  double computeParam(Id ivar, Id jvar) const;

  MatrixDense simulateSpectralOmega(Id nb) const override;

  double evalSpectrum(const VectorDouble& freq,
                      Id ivar = 0,
                      Id jvar = 0) const override;

  double evalSpectrumRatio(const VectorDouble& freq,
                           Id ivar,
                           Id jvar,
                           const ACov* cov0 = nullptr) const override;

protected:
  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;
  void _optimizationSetTarget(SpacePoint& pt) const override;

private:
  void _optimizationPreProcess(Id mode, const std::vector<SpacePoint>& ps) const override;
  void _optimizationPostProcess() const override;

private:
  Id _nVar;
  std::shared_ptr<const CorAniso> _corRef;     // the reference covariance function (ivar = 0)
  mutable CorAniso _corMatern; // an auxiliary Matern covariance function
  VectorDouble _coeffScales;   // scale factor for each variable
  VectorDouble _params;        // Matern parameter for each variable
  MatrixSymmetric _C0;         // the covariance matrix
  MatrixSymmetric _Rr;         // the scale factor matrix
  MatrixSymmetric _Nu;         // the parameter parameter matrix
  MatrixSymmetric _Tau;        // the maximal covariance matrix (_C0/_Tau should be positive definite)
};
} // namespace gstlrn
