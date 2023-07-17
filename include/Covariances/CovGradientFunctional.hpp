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
#include "Covariances/ACovGradient.hpp"
#include "Covariances/CovContext.hpp"

class Rotation;

/**
 * Class dedicated to manipulating a variables and its derivatives.
 * This feature is limited to the monovariate case
 */
class GSTLEARN_EXPORT CovGradientFunctional: public ACovGradient
{
public:
  CovGradientFunctional(const ECov& type, const CovContext& ctxt);
  CovGradientFunctional(const CovGradientFunctional& r);
  CovGradientFunctional(const CovAniso& r);
  CovGradientFunctional& operator=(const CovGradientFunctional& r);
  virtual ~CovGradientFunctional();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovGradientFunctional)

  void evalZAndGradients(const SpacePoint& p1,
                         const SpacePoint& p2,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode* mode = nullptr,
                         bool flagGrad = false) const override;

private:
  void   _calculateTrTtr(const VectorDouble& d,
                         VectorDouble& u,
                         VectorDouble& trttr) const;
};

