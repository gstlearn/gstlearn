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

#include "Basic/Vector.hpp"

#include "Basic/Tensor.hpp"
#include "Basic/AException.hpp"
#include "Basic/IClonable.hpp"
#include "Covariances/ACovGradient.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/ECov.hpp"

class Rotation;

/**
 * Class dedicated to manipulating a variables and its derivatives.
 * This feature is limited to the monovariate case
 */
class CovGradientFunctional: public ACovGradient
{
public:
  CovGradientFunctional(const ECov& type, const CovContext& ctxt);
  CovGradientFunctional(const CovGradientFunctional& r);
  CovGradientFunctional& operator=(const CovGradientFunctional& r);
  virtual ~CovGradientFunctional();

  virtual IClonable* clone() const override;

  void evalZAndGradients(const SpacePoint& p1,
                         const SpacePoint& p2,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGg,
                         const CovCalcMode& mode = CovCalcMode(),
                         bool flagGrad = false) const override;

private:
  void   _calculateTrTtr(const VectorDouble& d,
                         VectorDouble& u,
                         VectorDouble& trttr) const;
};

