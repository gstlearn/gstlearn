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
#include "LinearOp/InvNuggetOp.hpp"
#include "LinearOp/PrecisionOpMultiMatrix.hpp"
#include "LinearOp/ProjMultiMatrix.hpp"
#include "Matrix/MatrixSparse.hpp"
#include <memory>

namespace gstlrn
{
SPDEOpMatrix::SPDEOpMatrix(const PrecisionOpMultiMatrix* pop,
                           const ProjMultiMatrix* A,
                           const InvNuggetOp* invNoise)
  : SPDEOp(pop, A, invNoise, nullptr, nullptr)
  , _QpAinvNoiseAt(std::make_shared<MatrixSparse>(0, 0))
  , _chol(nullptr)
{
  _QpAinvNoiseAt->resize(pop->getSize(), pop->getSize());
  if (A != nullptr)
  {
    _QpAinvNoiseAt->prodNormMatMatInPlace(A->getProj(), invNoise, true);
  }
  _QpAinvNoiseAt->addMat(*pop->getQ());
}

SPDEOpMatrix::~SPDEOpMatrix()
{
  delete _chol;
}

Id SPDEOpMatrix::_solve(const constvect inv, vect outv) const
{
  if (_chol == nullptr)
    _chol = new CholeskySparse(*_QpAinvNoiseAt);
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
Id SPDEOpMatrix::_addToDest(const constvect inv, vect outv) const
{
  return _QpAinvNoiseAt->addToDest(inv, outv);
}

double SPDEOpMatrix::computeLogDetOp(Id nbsimu) const
{
  DECLARE_UNUSED(nbsimu);

  if (_chol == nullptr)
    _chol = new CholeskySparse(*_QpAinvNoiseAt); // TODO avoid to do it twice
  return _chol->computeLogDeterminant();
}

/**
 * @brief Computing Standard deviation of the estimation error
 * using partial_invert of a Sparse Cholesky matrix
 *
 * @param dat Vector of Data
 * @param nMC  Number of Monte-Carlo simulations (unused)
 * @param seed Random seed for the Monte-Carlo simulations (unused)
 * @param projK Projection Matrix used for Kriging
 * @param projS Projection matrix used for Simulations (unused)
 * @return VectorDouble
 */
VectorDouble SPDEOpMatrix::stdev(const VectorDouble& dat,
                                 Id nMC,
                                 Id seed,
                                 const ProjMulti* projK,
                                 const ProjMulti* projS) const
{
  DECLARE_UNUSED(dat);
  DECLARE_UNUSED(nMC);
  DECLARE_UNUSED(seed);
  DECLARE_UNUSED(projS);

  if (_chol == nullptr)
    _chol = new CholeskySparse(*_QpAinvNoiseAt); // TODO avoid to do it twice

  const auto* proj            = dynamic_cast<const ProjMultiMatrix*>(projK);
  const MatrixSparse* projmat = proj->getProj();

  VectorDouble result(projmat->getNRows());
  _chol->stdev(result, projmat, true); // true for standard deviation

  return result;
}
} // namespace gstlrn