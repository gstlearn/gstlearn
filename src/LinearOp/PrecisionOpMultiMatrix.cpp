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
{
  buildQop();
  _prepareMatrix();
}

MatrixSparse PrecisionOpMultiMatrix::_prepareMatrixStationary(int icov, const MatrixSparse* Q)
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
      MatrixSparse::glueInPlace(&currentCol,&copy ,1,0);
    }
    MatrixSparse::glueInPlace(&current, &currentCol, 0, 1);
  }
  return current;
}

MatrixSparse PrecisionOpMultiMatrix::_prepareMatrixNoStat(int icov, const MatrixSparse* Q)
{
  int n = PrecisionOpMulti::size(icov);
  MatrixSparse diag(n,n);
  const MatrixSparse empty(n,n);

  MatrixSparse bigQ = MatrixSparse(0,0);
  for (int jvar = 0; jvar < _getNVar(); jvar++)
  {
    MatrixSparse::glueInPlace(&bigQ, Q, 1, 1);
  }

  MatrixSparse bigLambda = MatrixSparse(0,0);
  int s = 0;
  for (int jvar = 0; jvar < _getNVar(); jvar++)
  {
    MatrixSparse currentCol(0,0);
    for (int ivar = 0; ivar< _getNVar(); ivar++)
    {
      if (ivar < jvar)
      {
        MatrixSparse::glueInPlace(&currentCol,&empty ,1,0);
      }
      else 
      {
        diag.setDiagonal(_invCholSillsNoStat[icov][s++]);
        MatrixSparse::glueInPlace(&currentCol,&diag ,1,0);
      }
    }
    MatrixSparse::glueInPlace(&bigLambda, &currentCol, 0, 1);
  }

  MatrixSparse result;
  result.prodNormMatMatInPlace(bigLambda,bigQ,true);
  return result;
}

void PrecisionOpMultiMatrix::_prepareMatrix()
{
  MatrixSparse current;
  for (int istruct = 0; istruct < _getNCov(); istruct++)
  {   
    const MatrixSparse* Q = ((PrecisionOpCs*)_pops[istruct])->getQ();

    if (_model->getVariableNumber() == 1)
    {
      MatrixSparse::glueInPlace(this,Q,1,1);
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
      MatrixSparse::glueInPlace(this, &current,1,1);
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

int PrecisionOpMultiMatrix::_addToDest(const Eigen::VectorXd& inv,
                                      Eigen::VectorXd& outv) const
{
  return MatrixSparse::_addToDest(inv,outv);
}

                                        