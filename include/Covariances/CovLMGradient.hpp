/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
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

  void evalZAndGradients(const SpacePoint& p1,
                         const SpacePoint& p2,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode& mode = CovCalcMode(),
                         bool flagGrad = false) const;
  void evalZAndGradients(const VectorDouble& vec,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode& mode = CovCalcMode(),
                         bool flagGrad = false) const;

private:
  void _initGradients(double& covVal,
                      VectorDouble& covGp,
                      VectorDouble& covGG,
                      bool flagGrad) const;
};
