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
#include "Basic/AStringable.hpp"
#include "Covariances/CovAniso.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSymmetric.hpp"

namespace gstlrn{

PrecisionOpMultiMatrix::PrecisionOpMultiMatrix(Model* model,
                                   const VectorMeshes& meshes)
  : PrecisionOpMulti(model,meshes,false,false)
  , _Q(MatrixSparse(0,0))
{
  buildQop(false);
  _prepareMatrix();
}

MatrixSparse PrecisionOpMultiMatrix::_prepareMatrixStationary(Id icov, const MatrixSparse* Q) const
{
  MatrixSymmetric sills = _sills[icov];
  sills.invert();
    
  MatrixSparse current(0,0);
  for (Id jvar = 0; jvar < _getNVar(); jvar++)
  {
    MatrixSparse currentCol(0,0);
    for (Id ivar = 0; ivar < _getNVar(); ivar++)
    {
      MatrixSparse copy = *Q;
      copy.prodScalar(sills.getValue(ivar,jvar));
      MatrixSparse::glueInPlace(&currentCol,&copy ,1,0);
    }

    MatrixSparse::glueInPlace(&current, &currentCol, 0, 1);
  }
  return current;
}

MatrixSparse PrecisionOpMultiMatrix::_prepareMatrixNoStat(Id icov, const MatrixSparse* Q) const
{
  Id n = PrecisionOpMulti::size(icov);
  auto nvar = _getNVar();
  const MatrixSparse empty(n,n);
  MatrixSparse diag(n,n);

  MatrixSparse bigQ(0,0);
  for (Id jvar = 0; jvar < nvar; jvar++)
  {
    MatrixSparse::glueInPlace(&bigQ, Q, 1,1);
  }
  MatrixSparse bigLambda(0,0);
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    MatrixSparse currentRow(0,0);
    for (Id jvar = 0; jvar < nvar; jvar++)
    {
      if (jvar <= ivar)
      {
        const auto& vec = &_invCholSillsNoStat[icov][IND(ivar, jvar, nvar)];
        constvect vecs(vec->data(),vec->size());
        diag.setDiagonal(vecs);
        MatrixSparse::glueInPlace(&currentRow,&diag ,1,0);
      }
      else 
      {
        MatrixSparse::glueInPlace(&currentRow,&empty ,1,0);
      }
    }
    MatrixSparse::glueInPlace(&bigLambda, &currentRow, 0, 1);
  }

  MatrixSparse result(bigQ.getNRows(),bigQ.getNCols());
  result.prodNormMatMatInPlace(&bigLambda,&bigQ,false);
  return result;
}

const MatrixSparse* PrecisionOpMultiMatrix::getQ() const
{
  if (_isSingle())
  {
    return ((PrecisionOpMatrix*)_pops[0])->getQ();
  }
  return &_Q;
}

void PrecisionOpMultiMatrix::_prepareMatrix()
{
  if (_isSingle()) return;

  MatrixSparse current(0, 0);
  for (Id istruct = 0; istruct < _getNCov(); istruct++)
  {
    const MatrixSparse* Q = ((PrecisionOpMatrix*)_pops[istruct])->getQ();

    if (_model->getNVar() == 1)
    {
      MatrixSparse::glueInPlace(&_Q, Q, 1, 1);
    }
    else
    {
      if (_isNoStatForVariance[istruct])
      {
        current = _prepareMatrixNoStat(istruct, Q);
      }
      else
      {
        current = _prepareMatrixStationary(istruct, Q);
      }
      MatrixSparse::glueInPlace(&_Q, &current, 1, 1);
    }
  }
}

PrecisionOpMultiMatrix::~PrecisionOpMultiMatrix()
{

}

void PrecisionOpMultiMatrix::_buildQop(bool stencil)
{
  if (stencil)
  {
    messerr("PrecisionOpMultiMatrix does not support stencil option\n");
  }
  for (Id icov = 0, number = _getNCov(); icov < number; icov++)
  {
    CovAniso* cova = _model->getCovAniso(_getCovInd(icov));
    _pops.push_back(new PrecisionOpMatrix(_meshes[icov], cova));
  }
}

Id PrecisionOpMultiMatrix::_addToDest(const constvect vecin, vect vecout) const
{
  return getQ()->addToDest(vecin, vecout);
}
}