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
#include "Covariances/ACov.hpp"
#include "Covariances/ECov.hpp"
#include "Covariances/ACovFunc.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovContext.hpp"

class Rotation;

/**
 * Class dedicated to manipulating a variables and its derivatives.
 * This feature is limited to the monovariate case
 */
class ACovGradient: public CovAniso
{
public:
  ACovGradient(const ECov& type, const CovContext& ctxt);
  ACovGradient(const ACovGradient& r);
  ACovGradient& operator=(const ACovGradient& r);
  virtual ~ACovGradient();

  virtual void evalZAndGradients(const SpacePoint& p1,
                                 const SpacePoint& p2,
                                 double& covVal,
                                 VectorDouble& covGp,
                                 VectorDouble& covGg,
                                 const CovCalcMode& mode = CovCalcMode(),
                                 bool flagGrad = false) const = 0;

private:
  bool _validVariable(int ivar) const;
};
