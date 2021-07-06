/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "Basic/Vector.hpp"
#include "PrecisionOp.hpp"
#include "ProjMatrix.hpp"

class HessianOp : public ALinearOp {

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
  void _evalDirect(const VectorDouble& in, VectorDouble& out) const override;

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
