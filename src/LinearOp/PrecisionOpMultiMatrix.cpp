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

#include "LinearOp/PrecisionOpMultiMatrix.hpp"
#include "LinearOp/ALinearOp.hpp"
#include "Matrix/MatrixSparse.hpp"
#include <Eigen/src/Core/Matrix.h>



PrecisionOpMultiMatrix::PrecisionOpMultiMatrix(Model* model,
                                   const std::vector<const AMesh*>& meshes)
  : PrecisionOpMulti(model,meshes)
{
 
}

PrecisionOpMultiMatrix::~PrecisionOpMultiMatrix()
{
}


int PrecisionOpMultiMatrix::_addToDest(const Eigen::VectorXd& inv,
                                      Eigen::VectorXd& outv) const
{
  MatrixSparse::addToDest(inv, outv);
  return 0;
}

                                            
