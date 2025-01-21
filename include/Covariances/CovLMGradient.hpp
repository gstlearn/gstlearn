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

#include "Space/ASpace.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Covariances/CovAnisoList.hpp"

class ASpace;
class SpacePoint;
class CovGradientNumerical;
class CovCalcMode;

class GSTLEARN_EXPORT CovLMGradient : public CovAnisoList
{
public:
  CovLMGradient(const ASpaceSharedPtr& space = ASpaceSharedPtr());
  CovLMGradient(const CovLMGradient& r);
  CovLMGradient(const CovAnisoList& r);
  CovLMGradient& operator= (const CovLMGradient &r);
  virtual ~CovLMGradient();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovLMGradient)

  // Add an elementary covariance structure
  virtual void addCovAniso(const CovAniso* cov) override;

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
  void _loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral& mat,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode *mode = nullptr) const override
  {
    ACov::_loadAndAddEvalCovMatBiPointInPlace(mat, p1, p2, mode);
  }
  void addEval0CovMatBiPointInPlace(MatrixSquareGeneral& mat,
                                    const CovCalcMode* mode) const override
  {
    ACov::addEval0CovMatBiPointInPlace(mat, mode);
  }
  void _addEvalCovMatBiPointInPlace(MatrixSquareGeneral& mat,
                                    const SpacePoint& pwork1,
                                    const SpacePoint& pwork2,
                                    const CovCalcMode* mode = nullptr) const override
  {
    ACov::_addEvalCovMatBiPointInPlace(mat, pwork1, pwork2, mode);
  }
  void _optimizationSetTarget(const SpacePoint& pt) const override
  {
    ACov::_optimizationSetTarget(pt); // TODO: cannot replace by CovAnisoList???
  }

private:
  static void _initGradients(double& covVal,
                             VectorDouble& covGp,
                             VectorDouble& covGG,
                             bool flagGrad);
};
