/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/ACovAnisoList.hpp"

class ASpace;
class SpacePoint;
class CovGradientNumerical;
class CovCalcMode;

class GSTLEARN_EXPORT CovLMGradient : public ACovAnisoList
{
public:
  CovLMGradient(const ASpace* space = nullptr);
  CovLMGradient(const CovLMGradient &r);
  CovLMGradient& operator= (const CovLMGradient &r);
  virtual ~CovLMGradient();

  virtual IClonable* clone() const override { return new CovLMGradient(*this); };

  // Add an elementary covariance structure
  virtual void addCov(const CovAniso* cov) override;

  void evalZAndGradients(const SpacePoint& p1,
                         const SpacePoint& p2,
                         double* covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGg,
                         const CovCalcMode& mode = CovCalcMode(),
                         bool flagGrad = false) const;
  void evalZAndGradients(const VectorDouble& vec,
                         double* covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGg,
                         const CovCalcMode& mode = CovCalcMode(),
                         bool flagGrad = false) const;

private:
  void _initGradients(double* covVal,
                      VectorDouble& covGp,
                      VectorDouble& covGg,
                      bool flagGrad) const;
};
