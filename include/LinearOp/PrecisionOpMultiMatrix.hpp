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

#include "LinearOp/ALinearOp.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "gstlearn_export.hpp"
#include "Model/Model.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "Basic/VectorNumT.hpp"

#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
#endif
class Model;

/**
 * Class to store objects for SPDE
 */
class GSTLEARN_EXPORT PrecisionOpMultiMatrix :  public PrecisionOpMulti, public MatrixSparse
{
public:
  PrecisionOpMultiMatrix(Model* model = nullptr, 
                   const std::vector<const AMesh*>& meshes = std::vector<const AMesh*>());
  PrecisionOpMultiMatrix(const PrecisionOpMulti &m)= delete;
  PrecisionOpMultiMatrix& operator= (const PrecisionOpMulti &m)= delete;
  virtual ~PrecisionOpMultiMatrix();

  #ifndef SWIG
  
  protected:
    int    _addToDest(const Eigen::VectorXd& inv,
                          Eigen::VectorXd& outv) const override;

  #endif

};
