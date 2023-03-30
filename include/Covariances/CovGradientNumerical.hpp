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

#include "Basic/Tensor.hpp"
#include "Basic/AException.hpp"
#include "Basic/ICloneable.hpp"
#include "Covariances/ACovGradient.hpp"
#include "Covariances/CovContext.hpp"

class Rotation;

/**
 * Class dedicated to manipulating a variables and its derivatives.
 * This feature is limited to the monovariate case
 */
class GSTLEARN_EXPORT CovGradientNumerical: public ACovGradient
{
public:
  CovGradientNumerical(const ECov& type, double ballRadius, const CovContext& ctxt);
  CovGradientNumerical(const CovGradientNumerical& r);
  CovGradientNumerical& operator=(const CovGradientNumerical& r);
  virtual ~CovGradientNumerical();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovGradientNumerical)

  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode& mode = CovCalcMode()) const override;
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const override;

  virtual double getBallRadius() const override { return _ballRadius; }

  void evalZAndGradients(const SpacePoint& p1,
                         const SpacePoint& p2,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode& mode = CovCalcMode(),
                         bool flagGrad = false) const override;

private:
  double _evalZZ(int ivar,
                 int jvar,
                 const SpacePoint& p1,
                 const SpacePoint& p2,
                 const CovCalcMode& mode = CovCalcMode()) const;
  double _evalZGrad(int ivar,
                    int jvar,
                    int idim,
                    const SpacePoint& p1,
                    const SpacePoint& p2,
                    const CovCalcMode& mode = CovCalcMode()) const;
  double _evalGradGrad(int ivar,
                       int jvar,
                       int idim,
                       int jdim,
                       const SpacePoint& p1,
                       const SpacePoint& p2,
                       const CovCalcMode& mode = CovCalcMode()) const;

private:
  double _ballRadius;   /*! Radius of the Ball for Numerical Gradient calculation */
};

