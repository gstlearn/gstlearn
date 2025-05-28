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
#include "LinearOp/MatrixSymmetricSim.hpp"
#include "Matrix/MatrixSparse.hpp"

SPDEOpMatrix::SPDEOpMatrix(const PrecisionOpMultiMatrix* pop,
                           const ProjMultiMatrix* A,
                           const MatrixSparse* invNoise,
                           const ProjMultiMatrix* projOut)
  : SPDEOp(pop,
           A,
           (invNoise == nullptr) ? nullptr : new MatrixSymmetricSim(invNoise),
           nullptr,
           nullptr,
           projOut, 
           projOut,
           true)
  , _QpAinvNoiseAt(MatrixSparse(0, 0))
  , _chol(nullptr)
{
  _QpAinvNoiseAt.resize(pop->getSize(), pop->getSize());
  if (A != nullptr)
  {
    _QpAinvNoiseAt.prodNormMatMatInPlace(A->getProj(), invNoise, true);
  }
  _QpAinvNoiseAt.addMatInPlace(*pop->getQ());
}

SPDEOpMatrix::~SPDEOpMatrix()
{
  delete _chol; 
}

int SPDEOpMatrix::_solve(const constvect inv, vect outv) const
{
  if (_chol == nullptr)
    _chol = new CholeskySparse(&_QpAinvNoiseAt);
  return _chol->solve(inv, outv);
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
int SPDEOpMatrix::_addToDest(const constvect inv, vect outv) const
{
 return _QpAinvNoiseAt.addToDest(inv,outv);
}

double SPDEOpMatrix::computeLogDetOp(int nbsimu) const
{
  DECLARE_UNUSED(nbsimu);

  if (_chol == nullptr)
    _chol = new CholeskySparse(&_QpAinvNoiseAt);
  return _chol->computeLogDeterminant();
}
