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
#include "LinearOp/PrecisionOpMultiMatrix.hpp"
#include "LinearOp/ProjMultiMatrix.hpp"

SPDEOpMatrix::SPDEOpMatrix(const PrecisionOpMultiMatrix* pop, const ProjMultiMatrix* A, const MatrixSparse* invNoise)
: SPDEOp(pop,A,invNoise)
{
  _QpAinvNoiseAt = invNoise->clone();
}

SPDEOpMatrix::~SPDEOpMatrix() {}

int SPDEOpMatrix::getSize() const
{ 
  return _Q->getSize(); 
}
/*****************************************************************************/
/*!
**  Evaluate the product (by the SPDEOpMatrix) : 'outv' = I * 'inv' = 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
int SPDEOpMatrix::_addToDest(const Eigen::VectorXd& inv,
                          Eigen::VectorXd& outv) const
{
  for (int i = 0, n = getSize(); i < n; i++)
    outv[i] += inv[i];
  return 0;
}
