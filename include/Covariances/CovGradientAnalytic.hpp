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

#include "Covariances/CovGradientGeneric.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class ACov;
class CovAniso;
/**
 * \brief
 * This class describes the Covariance to be used when processing
 * Data, Gradients and Tangents in the Potential Method framework.
 * This covariance is based on the initial covariance of the Data and derives
 * all simple and cross covariances for using Gradients and/or Tangents.
 * It uses Functional Derivation and therefore is suitable for a limited
 * set of differentiable covariances.
 */
class GSTLEARN_EXPORT CovGradientAnalytic: public CovGradientGeneric
{
public:
  CovGradientAnalytic(const CovAniso& cova);
  CovGradientAnalytic(const CovGradientAnalytic& r);
  CovGradientAnalytic& operator=(const CovGradientAnalytic& r) = delete;
  virtual ~CovGradientAnalytic();

  /// ICloneable Interface
  IMPLEMENT_CLONING(CovGradientAnalytic)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;
  void setFlagCalculateGG(bool status) { _flagCalculateGG = status; }

protected:
  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;
  void _optimizationSetTarget(SpacePoint& pt) const override;

private:
  bool _isValid() const override;
  void _calculateTrTtr() const;
  void _evalZAndGradients(const SpacePoint& p1, const SpacePoint& p2) const;
  const CovAniso* _getCovRefAniso() const;
  bool _onePointHasChanged(const SpacePoint& p1, const SpacePoint& p2) const;

private:
  mutable bool _flagCalculateGG;

  // Memory SpacePoints... used to check if bipoint has changed
  mutable SpacePoint _p1;
  mutable SpacePoint _p2;

  // covpp  Covariance value
  mutable double _covpp;
  // covGp  Covariance <G[i](x0+x,y0+y,z0+z), P(x0,y0,z0)> (dim=3)
  mutable VectorDouble _covGp;
  // covGG  Covariance <G[i](x0+x,y0+y,z0+z), G[j](x0,y0,z0)> (dim=3)
  mutable VectorDouble _covGG;

  // Working arrays
  mutable VectorDouble _dF;
  mutable VectorDouble _uF;
  mutable VectorDouble _hF;
  mutable VectorDouble _Tr;
  mutable VectorDouble _trttr;
};

} // namespace gstlrn
