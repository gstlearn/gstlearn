/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/ECov.hpp"

#include "Basic/Tensor.hpp"
#include "Basic/AException.hpp"
#include "Basic/ICloneable.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/ACovFunc.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovContext.hpp"

class Rotation;

/**
 * Class dedicated to manipulating a variables and its derivatives.
 * This feature is limited to the monovariate case
 */
class GSTLEARN_EXPORT ACovGradient: public CovAniso
{
public:
  ACovGradient(const ECov& type, const CovContext& ctxt);
  ACovGradient(const ACovGradient& r);
  ACovGradient(const CovAniso& r);
  ACovGradient& operator=(const ACovGradient& r);
  virtual ~ACovGradient();

  virtual void evalZAndGradients(const SpacePoint& p1,
                                 const SpacePoint& p2,
                                 double& covVal,
                                 VectorDouble& covGp,
                                 VectorDouble& covGG,
                                 const CovCalcMode& mode = CovCalcMode(),
                                 bool flagGrad = false) const = 0;
};
