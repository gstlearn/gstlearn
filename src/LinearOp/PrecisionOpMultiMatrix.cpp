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
#include "Covariances/CovAniso.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

PrecisionOpMultiMatrix::PrecisionOpMultiMatrix(Model* model,
                                   const VectorMeshes& meshes)
  : PrecisionOpMulti(model,meshes,false)
  , _Q(MatrixSparse(0,0))
{
  buildQop();
  _prepareMatrix();
}

MatrixSparse PrecisionOpMultiMatrix::_prepareMatrixStationary(int icov, const MatrixSparse* Q) const
{
  MatrixSquareSymmetric sills = *_invCholSillsStat[icov].getMatrix();
  sills.invert();
    
  MatrixSparse current = MatrixSparse(0,0);
  for (int jvar = 0; jvar < _getNVar(); jvar++)
  {
    MatrixSparse currentCol = MatrixSparse(0,0);
    for (int ivar = 0; ivar < _getNVar(); ivar++)
    {
      MatrixSparse copy = *Q;
      copy.prodScalar(sills.getValue(ivar,jvar));
      MatrixSparse::glueInPlace(&currentCol,&copy ,1,0);
    }

    MatrixSparse::glueInPlace(&current, &currentCol, 0, 1);
  }
  return current;
}

MatrixSparse PrecisionOpMultiMatrix::_prepareMatrixNoStat(int icov, const MatrixSparse* Q) const
{
  int n = PrecisionOpMulti::size(icov);
  int nvar = _getNVar();
  const MatrixSparse empty(n,n);
  MatrixSparse diag(n,n);

  MatrixSparse bigQ = MatrixSparse(0,0);
  for (int jvar = 0; jvar < nvar; jvar++)
  {
    MatrixSparse::glueInPlace(&bigQ, Q, 1,1);
  }
  MatrixSparse bigLambda = MatrixSparse(0,0);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    MatrixSparse currentRow(0,0);
    for (int jvar = 0; jvar < nvar; jvar++)
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
  for (int istruct = 0; istruct < _getNCov(); istruct++)
  {
    const MatrixSparse* Q = ((PrecisionOpMatrix*)_pops[istruct])->getQ();

    if (_model->getVariableNumber() == 1)
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

void PrecisionOpMultiMatrix::_buildQop()
{
  for (int icov = 0, number = _getNCov(); icov < number; icov++)
  {
    CovAniso* cova = _model->getCova(_getCovInd(icov));
    _pops.push_back(new PrecisionOpMatrix(_meshes[icov], cova));
  }
}

int PrecisionOpMultiMatrix::_addToDest(const constvect vecin,
                                           vect vecout) const
{
  return getQ()->addToDest(vecin, vecout);
}
