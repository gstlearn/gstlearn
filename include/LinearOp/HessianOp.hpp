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

#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
  #include <Eigen/src/Core/Matrix.h>
#endif

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


  void setLambda(const Eigen::VectorXd& lambda) 
  {
    for (int i = 0; i < (int)_lambda.size(); i++) 
      _lambda[i] = lambda[i]; 
  }

 void setLambda(const VectorDouble& lambda) 
  {
    for (int i = 0; i < (int)_lambda.size(); i++) 
      _lambda[i] = lambda[i]; 
  }

#ifndef SWIG
protected:
  int _addToDest(const constvect& inv,
                 vect& outv) const override;


private:
  bool                 _isInitialized;
  bool                 _flagSeismic;
  PrecisionOp*         _pMat;     // External pointer
  const ProjMatrix*    _projData; // External pointer
  const ProjMatrix*    _projSeis; // External pointer
  std::vector<double>         _indic;
  std::vector<double>          _propSeis;
  std::vector<double>          _varSeis;
  std::vector<double>          _lambda;
  mutable std::vector<double>  _workp;
  mutable std::vector<double>  _workx;
  mutable std::vector<double>  _workv;
  mutable std::vector<double>  _works;
#endif
};

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(HessianOp)
#endif