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

#include "LinearOp/MatrixSquareSymmetricSim.hpp"
#include "Basic/AStringable.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "LinearOp/ACholesky.hpp"
#include <Eigen/src/Core/Matrix.h>

MatrixSquareSymmetricSim::MatrixSquareSymmetricSim(const AMatrix* m,
                                                   bool inverse)
  : ASimulable()
  , _inverse(inverse)
  , _factor(nullptr)

{
  if (!m->isSquare())
  {
    messerr("The matrix must be square!");
    return;
  }

  if (m->isSparse())
  {
    const MatrixSparse* matCS = dynamic_cast<const MatrixSparse*>(m);
    if (matCS != nullptr) _factor = new CholeskySparse(matCS);
  }
  else
  {
    const MatrixSquareSymmetric* matSym = dynamic_cast<const MatrixSquareSymmetric*>(m);
    if (matSym != nullptr) _factor = new CholeskyDense(matSym);
  }
  if (_factor == nullptr)
  {
    messerr("The Input matrix is not valid");
    messerr("It should be either:");
    messerr("- a MatrixSparse");
    messerr("- a MatrixSquareSymmetric");
    return;
  }
}

MatrixSquareSymmetricSim::~MatrixSquareSymmetricSim()
{
  delete _factor;
  _factor = nullptr;
}

int MatrixSquareSymmetricSim::_addToDest(const constvect inv, vect outv) const
{
  if (_inverse) return _factor->getMatrix()->addProdMatVecInPlace(inv, outv);
  return _factor->addSolveX(inv, outv);
}

int MatrixSquareSymmetricSim::_addSimulateToDest(const constvect whitenoise,
                                                 vect outv) const
{
  if (_inverse) return _factor->addInvLtX(whitenoise, outv);
  return _factor->addLX(whitenoise, outv);
}

const AMatrix* MatrixSquareSymmetricSim::getMatrix() const
{
  if (_factor == nullptr) return nullptr;
  return _factor->getMatrix();
}

int MatrixSquareSymmetricSim::getSize() const
{
  if (_factor == nullptr) return 0;
  return _factor->getSize();
}