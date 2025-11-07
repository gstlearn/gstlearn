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
  CorGneiting(const ECov& type, const CovContext& ctxt);
  CorGneiting(const CorAniso* covS, const CorAniso* covT, double separability = 1.0);
  CorGneiting(const CorGneiting& r);
  CorGneiting& operator=(const CorGneiting& r);
  virtual ~CorGneiting();
  static CorGneiting* create(const ECov& type,
                             const CovContext& ctxt,
                             const VectorDouble& params,
                             const VectorDouble& ranges = VectorDouble(),
                             const VectorDouble& angles = VectorDouble(),
                             double separability        = 1.0,
                             bool flagRange             = false);
  IMPLEMENT_CLONING(CorGneiting)

  bool isConsistent(const ASpace* space) const override
  {
    DECLARE_UNUSED(space)
    return true;
  }

  /// ACov Interface

  Id getNVar() const override { return 1; }
  MatrixDense simulateSpectralOmega(Id nb) const override;

protected:
  // void _optimizationSetTarget(SpacePoint& pt) const override;
  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;
  // private:
  //  void _optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const override;
  // void _optimizationPostProcess() const override;

private:
  std::shared_ptr<const CorAniso> _corS;
  std::shared_ptr<const CorAniso> _corT;
  double _separability;
  mutable CorAniso _corSCopy;
};

} // namespace gstlrn