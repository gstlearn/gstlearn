/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "PrecisionOp.hpp"
#include "ProjMatrix.hpp"

class GSTLEARN_EXPORT HessianOp : public ALinearOp {

public:
	HessianOp();
	virtual ~HessianOp();

  int  init(PrecisionOp*  pmat,
            const ProjMatrix*   projdata,
            const ProjMatrix*   projseis,
            const VectorDouble& indic,
            const VectorDouble& propseis,
            const VectorDouble& varseis);
  /*!  Returns the dimension of the matrix */
  int  getSize() const override { return _pMat->getSize(); }
  /*!  Set the initial vector */
  void setLambda(const VectorDouble& lambda) { _lambda = lambda; };

protected:
  void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;

private:
  bool                 _isInitialized;
  bool                 _flagSeismic;
  PrecisionOp*   _pMat;
  const ProjMatrix*    _projData;
  const ProjMatrix*    _projSeis;
  VectorDouble         _indic;
  VectorDouble         _propSeis;
  VectorDouble         _varSeis;
  VectorDouble         _lambda;
  mutable VectorDouble _workp;
  mutable VectorDouble _workx;
  mutable VectorDouble _workv;
  mutable VectorDouble _works;
};
