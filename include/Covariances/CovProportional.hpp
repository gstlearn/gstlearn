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
#include "Covariances/CovBase.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
class ACor;
class AFunctional;
class CovInternal;

class GSTLEARN_EXPORT CovProportional: public CovBase
{
public:
    
  CovProportional(ACor* cor = nullptr,const MatrixSquareSymmetric &sills = MatrixSquareSymmetric());
  CovProportional(const CovProportional &r) = delete;
  CovProportional& operator=(const CovProportional &r) = delete;
  virtual ~CovProportional();

  void setCor(ACor* cor) override;
protected:
    
  virtual void _addEvalCovMatBiPointInPlace(MatrixSquareGeneral& mat,
                                            const SpacePoint& p1,
                                            const SpacePoint& p2,
                                            const CovCalcMode* mode = nullptr) const override;
 
protected:
    mutable MatrixSquareGeneral _workMat;
};
