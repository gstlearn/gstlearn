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

#include "Covariances/CovContext.hpp"
#include "Space/ASpace.hpp"
#include "gstlearn_export.hpp"

#include "Covariances/CovAnisoList.hpp"

class ASpace;
class SpacePoint;
class CovGradientNumerical;
class CovCalcMode;

class GSTLEARN_EXPORT CovLMGradient : public CovAnisoList
{
public:
  CovLMGradient(const CovContext& ctxt = CovContext());
  CovLMGradient(const CovLMGradient& r);
  CovLMGradient(const CovAnisoList& r);
  CovLMGradient& operator= (const CovLMGradient &r);
  virtual ~CovLMGradient();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovLMGradient)

  // Add an elementary covariance structure
  virtual void addCov(const CovBase* cov) override;

  /// ACov interface
  
  void evalZAndGradients(const SpacePoint& p1,
                         const SpacePoint& p2,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode* mode = nullptr,
                         bool flagGrad = false) const;
  void evalZAndGradients(const VectorDouble& vec,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode* mode = nullptr,
                         bool flagGrad = false) const;

protected:
  void _optimizationSetTarget(SpacePoint& pt) const override
  {
    ACov::_optimizationSetTarget(pt); 
  }

private:
  static void _initGradients(double& covVal,
                             VectorDouble& covGp,
                             VectorDouble& covGG,
                             bool flagGrad);
};
