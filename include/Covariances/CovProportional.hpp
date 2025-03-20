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
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
class ACor;
class AFunctional;
class CovInternal;

class GSTLEARN_EXPORT CovProportional: public CovBase
{
public:
  CovProportional(ACov* cor = nullptr,const MatrixSquareSymmetric &sills = MatrixSquareSymmetric());
  CovProportional(const CovProportional &r);
  CovProportional& operator=(const CovProportional &r);
  virtual ~CovProportional();

  void setCor(ACov* cor) override;
  IMPLEMENT_CLONING(CovProportional)
protected:
  double _eval(const SpacePoint& p1, 
               const SpacePoint& p2,
               int ivar = 0, 
               int jvar = 0, 
               const CovCalcMode* mode = nullptr) const override;
protected:
    mutable MatrixSquareGeneral _workMat;
};
