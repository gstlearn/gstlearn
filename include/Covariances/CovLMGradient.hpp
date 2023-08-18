/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Covariances/ACovAnisoList.hpp"

class ASpace;
class SpacePoint;
class CovGradientNumerical;
class CovCalcMode;

class GSTLEARN_EXPORT CovLMGradient : public ACovAnisoList
{
public:
  CovLMGradient(const ASpace* space = nullptr);
  CovLMGradient(const CovLMGradient& r);
  CovLMGradient(const ACovAnisoList& r);
  CovLMGradient& operator= (const CovLMGradient &r);
  virtual ~CovLMGradient();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovLMGradient)

  // Add an elementary covariance structure
  virtual void addCov(const CovAniso* cov) override;

  /// ACov interface
  virtual void eval0MatInPlace(MatrixSquareGeneral &mat,
                               const CovCalcMode *mode = nullptr) const override;
  virtual void evalMatInPlace(const SpacePoint &p1,
                              const SpacePoint &p2,
                              MatrixSquareGeneral &mat,
                              const CovCalcMode *mode = nullptr) const override;
  /// Tell if the use of Optimization is enabled or not
  virtual bool isOptimEnabled() const override { return false; }

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


private:
  void _initGradients(double& covVal,
                      VectorDouble& covGp,
                      VectorDouble& covGG,
                      bool flagGrad) const;
};
