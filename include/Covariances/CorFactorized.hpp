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

#include "Covariances/ACov.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Model/AModelFitSills.hpp"
#include "Enum/ESimuType.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <vector>

namespace gstlrn
{
  class ACov;
  class CovContext;

  /**
   * \brief
   * This class describes the **Covariance** as a list of elementary covariances defined on sub spaces
   * where the calculation rule is simple: the returned value is the **product** of each elementary covariance function
   * evaluated between the projection of the two points on the attached subspace.
   */

  class GSTLEARN_EXPORT CorFactorized: public ACov
  {
  public:
    CorFactorized(const CovContext& ctxt = CovContext());
    CorFactorized(const CorFactorized& r);
    CorFactorized& operator=(const CorFactorized& r);
    virtual ~CorFactorized();

    /// ICloneable Interface
    IMPLEMENT_CLONING(CorFactorized)

    static CorFactorized* create(const CovContext& ctxt);

    /// Interface for ACov
    double eval0(Id ivar = 0, Id jvar = 0, const CovCalcMode* mode = nullptr)
      const override;

    /// Interface for AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// CorFactorized Interface

    // Get the number of factors
    Id getNFac() const override;

    // Remove an elementary factor structure
    virtual void addFactor(const ACov& cov, const MatrixDense& proj);

    // Remove an elementary factor structure
    // void deleteFactor(Id ifac); Not used

    // Remove all elementary factor structures
    void deleteAllFactors();

    // Get the covariance structure of the ifac factor
    const ACov* getCovariance(Id ifac) const;

    // Get the projection matrix of the ifac factor
    MatrixDense getProjection(Id ifac = 0) const override;

    bool isValidForSimulation(const ESimuType& simuType) const override;
    SpectrumOnRN* simulateOnRN(Id ns = 1000) const override;

  protected:
    bool _isFactorIndexValid(Id ifac) const;
    double _eval(
      const SpacePoint& p1,
      const SpacePoint& p2,
      Id ivar = 0,
      Id jvar = 0,
      const CovCalcMode* mode = nullptr) const override;

#ifndef SWIG
    std::vector<std::shared_ptr<const ACov>>
      _covs; // the space covariance function (ndim, nvar >=1)
    std::vector<std::shared_ptr<const MatrixDense>>
      _projs; /// Vector of elementary projection to subspaces
#endif
  };
} // namespace gstlrn
