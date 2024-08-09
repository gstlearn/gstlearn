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

#include "gstlearn_export.hpp"

#include "LinearOp/ProjMatrix.hpp"

class PrecisionOp;

#ifndef SWIG
#  include "LinearOp/ALinearOpEigenCG.hpp"
DECLARE_EIGEN_TRAITS(HessianOp)
#else
#  include "LinearOp/ALinearOp.hpp"
#endif

class GSTLEARN_EXPORT HessianOp:
#ifndef SWIG
  public ALinearOpEigenCG<HessianOp>
#else
  public ALinearOp
#endif
{

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
  int  getSize() const override;
  /*!  Set the initial vector */
  void setLambda(const VectorDouble& lambda) { _lambda = lambda; };

#ifndef SWIG
protected:
  void _evalDirect(const Eigen::VectorXd& inv,
                   Eigen::VectorXd& outv) const override;

private:
  bool                 _isInitialized;
  bool                 _flagSeismic;
  PrecisionOp*         _pMat;     // External pointer
  const ProjMatrix*    _projData; // External pointer
  const ProjMatrix*    _projSeis; // External pointer
  VectorDouble         _indic;
  VectorDouble         _propSeis;
  VectorDouble         _varSeis;
  VectorDouble         _lambda;
  mutable VectorDouble _workp;
  mutable VectorDouble _workx;
  mutable VectorDouble _workv;
  mutable VectorDouble _works;
#endif
};

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(HessianOp)
#endif