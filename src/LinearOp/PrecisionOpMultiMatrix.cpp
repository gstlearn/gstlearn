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
#include <Eigen/src/Core/Matrix.h>
#include "LinearOp/PrecisionOpCs.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"



PrecisionOpMultiMatrix::PrecisionOpMultiMatrix(Model* model,
                                   const std::vector<const AMesh*>& meshes)
  : PrecisionOpMulti(model,meshes,false)
  , _Q(MatrixSparse(0,0))
{
  buildQop();
  _prepareMatrix();
}

MatrixSparse PrecisionOpMultiMatrix::_prepareMatrixStationary(int icov, const MatrixSparse* Q) const
{
  MatrixSquareSymmetric sills = _cholSills[icov];
  sills.invert();
    
  MatrixSparse current = MatrixSparse(0,0);
  for (int jvar = 0; jvar < _getNVar(); jvar++)
  {
    MatrixSparse currentCol = MatrixSparse(0,0);
    for (int ivar = 0; ivar < _getNVar(); ivar++)
    {
      MatrixSparse copy = *Q;
      copy.prodScalar(sills.getValue(ivar,jvar));
      int shift = ( ivar > 0 ) ? 1 : 0;
      MatrixSparse::glueInPlace(&currentCol,&copy ,shift,0);
    }
    int shift = ( jvar > 0 ) ? 1 : 0;

    MatrixSparse::glueInPlace(&current, &currentCol, 0, shift);
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
    int shift = ( jvar > 0 ) ? 1 : 0;
    MatrixSparse::glueInPlace(&bigQ, Q, shift, shift);
  }
  MatrixSparse bigLambda = MatrixSparse(0,0);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    MatrixSparse currentRow(0,0);
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      int shift = ( jvar > 0 ) ? 1 : 0;
      if (jvar <= ivar)
      {
        diag.setDiagonal(_invCholSillsNoStat[icov][IND(ivar,jvar,nvar)]);
        MatrixSparse::glueInPlace(&currentRow,&diag ,shift,0);
      }
      else 
      {
        MatrixSparse::glueInPlace(&currentRow,&empty ,shift,0);
      }
    }
    int shift = ( ivar > 0 ) ? 1 : 0;
    MatrixSparse::glueInPlace(&bigLambda, &currentRow, 0, shift);
  }

  MatrixSparse result(bigQ.getNRows(),bigQ.getNCols());
  result.prodNormMatMatInPlace(&bigLambda,&bigQ,false);
  return result;
}

const MatrixSparse* PrecisionOpMultiMatrix::getQ() const
{
  if (_isSingle())
  {
    return  ((PrecisionOpCs*)_pops[0])->getQ();
  }
    return &_Q;
  
}

void PrecisionOpMultiMatrix::_prepareMatrix()
{
  if (_isSingle()) return;

  MatrixSparse current(0,0);
  for (int istruct = 0; istruct < _getNCov(); istruct++)
  {   
    const MatrixSparse *Q = ((PrecisionOpCs*)_pops[istruct])->getQ();
    int shift = ( istruct > 0 ) ? 1 : 0;

    if (_model->getVariableNumber() == 1)
    {
      MatrixSparse::glueInPlace(&_Q,Q,shift,shift);
    }
    else 
    {
      if (_isNoStatForVariance[istruct])
      {
        current = _prepareMatrixNoStat(istruct,Q);
      }
      else 
      {
        current = _prepareMatrixStationary(istruct,Q);
      }
      MatrixSparse::glueInPlace(&_Q, &current,shift,shift);
    }
  }    
}

PrecisionOpMultiMatrix::~PrecisionOpMultiMatrix()
{

}

void PrecisionOpMultiMatrix::_buildQop()
{
  for (int i = 0, number = _getNCov(); i < number; i++)
  {
    _pops.push_back(new PrecisionOpCs(_meshes[i], _model, _getCovInd(i)));
  }
}

void PrecisionOpMultiMatrix::_makeReady()
{
  for (auto &e : _pops)
  {
    ((PrecisionOpCs*)e)->makeReady();
  }
}

int PrecisionOpMultiMatrix::_addToDestImpl(const Eigen::VectorXd &vecin,Eigen::VectorXd &vecout) const
{
  return getQ()->addToDest(vecin,vecout);
}
                                   