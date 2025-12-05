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
#include "Basic/ICloneable.hpp"
#include "Covariances/CovBase.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"

namespace gstlrn
{
class AFunctional;
class CovInternal;
class GSTLEARN_EXPORT CovProportional: public CovBase
{
public:
  CovProportional(ACov* cor = nullptr, const MatrixSymmetric& sills = MatrixSymmetric());
  CovProportional(const CovProportional& r);
  CovProportional& operator=(const CovProportional& r);
  virtual ~CovProportional();

  void setCor(ACov* cor) override;
  bool isValidForSpectral() const override;

  MatrixDense simulateSpectralOmega(Id ns) const override;
  IMPLEMENT_CLONING(CovProportional)

protected:
  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;

protected:
  mutable MatrixSquare _workMat;
};
} // namespace gstlrn