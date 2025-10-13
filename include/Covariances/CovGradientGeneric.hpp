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

#include "Covariances/ACov.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class ACov;
class CovAniso;
/**
 * \brief
 * This class describes the Covariance to be used when processing Data
 * and its gradient components.
 * This covariance is based on the initial covariance of the Data and derives
 * the simple and cross covariances of its gradient components.
 * It uses Numerical Derivation and therefore is suitable whatever the type
 * of covariance used for the Data variable.
 */
class GSTLEARN_EXPORT CovGradientGeneric: public ACov
{
public:
  CovGradientGeneric(const ACov& cova, double ballradius = 1);
  CovGradientGeneric(const CovGradientGeneric& r);
  CovGradientGeneric& operator=(const CovGradientGeneric& r) = delete;
  virtual ~CovGradientGeneric();

  /// ICloneable Interface
  IMPLEMENT_CLONING(CovGradientGeneric)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  double getBallRadius() const { return _ballRadius; }
  const ACov& getCovRef() const { return _covRef; }

protected:
  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;
  void _optimizationSetTarget(SpacePoint& pt) const override;

private:
  virtual bool _isValid() const;
  double _evalZGradientNumeric(const SpacePoint& p1,
                               const SpacePoint& p2,
                               Id idim,
                               const CovCalcMode* mode) const;
  double _evalGradientGradientNumeric(const SpacePoint& p1,
                                      const SpacePoint& p2,
                                      Id idim,
                                      Id jdim,
                                      const CovCalcMode* mode) const;

private:
  double _ballRadius;
  const ACov& _covRef; // Monovariate Reference structure
};

} // namespace gstlrn
