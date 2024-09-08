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
#include "LinearOp/SPDEOpMatrix.hpp"
#include "LinearOp/Cholesky.hpp"
#include "LinearOp/PrecisionOpMultiMatrix.hpp"
#include "LinearOp/ProjMultiMatrix.hpp"
#include "LinearOp/MatrixSquareSymmetricSim.hpp"
#include "Matrix/MatrixSparse.hpp"

SPDEOpMatrix::SPDEOpMatrix(const PrecisionOpMultiMatrix* pop, 
                           const ProjMultiMatrix* A, 
                           const MatrixSparse* invNoise)
: SPDEOp(pop,A,new MatrixSquareSymmetricSim(invNoise),1)
, _QpAinvNoiseAt(MatrixSparse(0,0))
{
  _QpAinvNoiseAt.resize(pop->getSize(), pop->getSize());
  _QpAinvNoiseAt.prodNormMatMatInPlace(A->getProj(),invNoise,true);
  _QpAinvNoiseAt.addMatInPlace(*pop->getQ());
}

SPDEOpMatrix::~SPDEOpMatrix()
{ 
}

int SPDEOpMatrix::_solve(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const
{
  _QpAinvNoiseAt.computeCholesky();
  return _QpAinvNoiseAt.solveCholesky(inv,outv);
}

/*****************************************************************************/
/*!
**  Evaluate the product (by the SPDEOpMatrix)
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
int SPDEOpMatrix::_addToDestImpl(const Eigen::VectorXd& inv,
                          Eigen::VectorXd& outv) const
{
 return _QpAinvNoiseAt.addToDest(inv,outv);
}
