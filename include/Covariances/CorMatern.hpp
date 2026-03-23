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
    CorMatern(const CovContext& ctxt,
              const ECov& type,
              const VectorDouble& params = VectorDouble(),
              const VectorDouble& kappas = VectorDouble(),
              const VectorDouble& ranges = VectorDouble(),
              const VectorDouble& angles = VectorDouble(),
              bool flagRange = true);
    CorMatern(const CorMatern& r);
    CorMatern& operator=(const CorMatern& r);
    virtual ~CorMatern();
    IMPLEMENT_CLONING(CorMatern)

    static CorMatern* create(const CovContext& ctxt,
                             const ECov& type,
                             const VectorDouble& params,
                             const VectorDouble& kappas,
                             const VectorDouble& ranges,
                             const VectorDouble& angles = VectorDouble(),
                             bool flagRange = true);

    const ECov& getType() const { return _corRef->getType(); }

    bool isConsistent(const ASpace* space) const override
    {
      return (space->getType() == ESpaceType::RN);
    }

    const CorAniso* getCorRef() const { return _corRef.get(); };

    /// ACov Interface
    Id getNVar() const override { return _nVar; }

    double getC0(Id ivar, Id jvar) const { return _C0.getValue(ivar, jvar); }

    double getNu(Id ivar, Id jvar) const { return _Nu.getValue(ivar, jvar); }

    double getKappa(Id ivar, Id jvar) const
    {
      return _Kappa.getValue(ivar, jvar);
    }

    static double computeParam(double nui, double nuj)
    {
      return 0.5 * (nui + nuj);
    }

    static double computeScale(double kappai, double kappaj)
    {
      return sqrt(0.5 * (kappai * kappai + kappaj * kappaj));
    }

    bool isValidForSimulation(const ESimuType& simuType) const override
    {
      return (getSpaceType() == ESpaceType::RN
              && simuType == ESimuType::SPECTRAL);
    }

    MatrixDense simulateSpectralOmega(Id nb) const override;
    SpectrumRN
      simulateSpectrumRN(Id ns, const ACov* cov0 = nullptr) const override;
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
                 Id ivar = 0,
                 Id jvar = 0,
                 const CovCalcMode* mode = nullptr) const override;
    void _optimizationSetTarget(SpacePoint& pt) const override;

  private:
    void
      _optimizationPreProcess(Id mode,
                              const std::vector<SpacePoint>& ps) const override;
    void _optimizationPostProcess() const override;

  private:
    Id _nVar;
    std::shared_ptr<const CorAniso>
      _corRef; // the reference covariance function (ivar = 0)
    mutable CorAniso _cor; // an auxiliary Matern covariance function
    MatrixSymmetric _Nu; // the parameter parameter matrix
    MatrixSymmetric _Kappa; // the scale factor matrix
    MatrixSymmetric _C0; // the covariance matrix
  };
} // namespace gstlrn
