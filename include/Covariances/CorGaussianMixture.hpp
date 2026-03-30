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
#include "Enum/ESpaceType.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include <memory>

namespace gstlrn
{
class CorAniso;
/**
 * \brief
 * This class describes the multivariate Matern correlation function.
 *
 */
class GSTLEARN_EXPORT CorGaussianMixture: public ACov
{
public:
  CorGaussianMixture(
    const CovContext& ctxt,
    const ECov& type,
    const VectorDouble& params = VectorDouble(),
    const VectorDouble& kappas = VectorDouble(),
    const VectorDouble& ranges = VectorDouble(),
    const VectorDouble& angles = VectorDouble(),
    bool flagRange             = true);
  CorGaussianMixture(const CorGaussianMixture& r);
  CorGaussianMixture& operator=(const CorGaussianMixture& r);
  virtual ~CorGaussianMixture();

  static CorGaussianMixture* create(
    const CovContext& ctxt,
    const ECov& type,
    const VectorDouble& params,
    const VectorDouble& kappas,
    const VectorDouble& ranges,
    const VectorDouble& angles = VectorDouble(),
    bool flagRange             = false);

  IMPLEMENT_CLONING(CorGaussianMixture)

  const ECov& getType() const { return _corRef->getType(); }

  bool isConsistent(const ASpace* space) const override
  {
    return (space->getType() == ESpaceType::RN);
  }

  const CorAniso* getCorRef() const { return _corRef.get(); }

  /// ACov Interface
  double getScaleDim(Id idim) const;
  VectorDouble getScales() const;
  void setScaleCor(Id idim, double scale);

  Id getNVar() const override { return _nVar; }
  double getC0(Id ivar, Id jvar) const { return _C0.getValue(ivar, jvar); }
  double getNu(Id ivar, Id jvar) const { return _Nu.getValue(ivar, jvar); }
  double getKappa(Id ivar, Id jvar) const { return _Kappa.getValue(ivar, jvar); }

  static double computeNu(double nui, double nuj) { return 0.5 * (nui + nuj); };
  double computeKappa(double kappai, double kappaj) const;

  // for Gneiting model...
  void setScaleGneiting(double val) { _scaleGneiting = val; }
  bool getScaleGneiting() const { return _scaleGneiting; }

  bool isValidForSimulation(const ESimuType& simuType) const override
  {
    return (getSpaceType() == ESpaceType::RN && simuType == ESimuType::SPECTRAL);
  }
  SpectrumOnRN* simulateOnRN(Id ns = 1000) const override;


protected:
  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;

private:
  Id _nVar;
  std::shared_ptr<const CorAniso> _corRef; // the reference covariance function (ivar = 0)
  mutable CorAniso _cor;                   // an auxiliary Matern or Cauchy covariance function
  MatrixSymmetric _Nu;                     // the parameter parameter matrix
  MatrixSymmetric _Kappa;                  // the scale factor matrix
  MatrixSymmetric _C0;                     // the covariance matrix
  double _scaleGneiting;                   // additional coefficient used for the implementation of the Gneiting model
};
} // namespace gstlrn
