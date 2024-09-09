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
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "LinearOp/Cholesky.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/AMatrixDense.hpp"

MatrixSquareSymmetricSim::MatrixSquareSymmetricSim(const AMatrix* m,bool inverse)
: _mat(m)
, _inverse(inverse)
, _empty(true)
, _factorDense(nullptr)
, _factorSparse(nullptr)
{
  if (m == nullptr) 
  { messerr("The matrix is null.");
    return;
  }
  if (!m->isSquare())
  {
    messerr("The matrix has to be square!");
    return;
  }

  _sparse = dynamic_cast<const MatrixSparse*>(m) !=nullptr;

  if (!_sparse)
  {
    bool symmetric = dynamic_cast<const MatrixSquareSymmetric*>(m) ==nullptr;
    if (!symmetric)
    {
          messerr("Warning, the matrix has to be symmetric.");
          messerr("You should create it as a MatrixSquareSymmetric class");
          messerr("instead of MatrixSquare or MatrixDense.");
    }
  }

  _empty = false;

}


MatrixSquareSymmetricSim::MatrixSquareSymmetricSim()
: _mat(nullptr) 
, _inverse(false)
, _empty(true)
, _factorDense(nullptr)
, _factorSparse(nullptr)
{
  _sparse = dynamic_cast<MatrixSparse*>(this) != nullptr;
}

int MatrixSquareSymmetricSim::_addToDest(const Eigen::VectorXd& inv,
                                               Eigen::VectorXd& outv) const
{  
  if (_inverse)
  {
    return _mat->addProdMatVecInPlace(inv,outv);
  }
   
  messerr("MatrixSquareSymmetricSim::_addToDest not implemented for inverse = false.");
  return 1;
}
int MatrixSquareSymmetricSim::_addSimulateToDest(const Eigen::VectorXd& whitenoise,
                                                       Eigen::VectorXd& outv) const
{  


  _prepare();
  if (isSparse())
  {
    if (isInverse())
         return _factorSparse->addSimulateToDest(whitenoise,outv);
    return _factorSparse->addToDest(whitenoise,outv);
  }
  
  if (!isSparse())
  {
    if (isInverse())
    {
      outv.noalias() += _factorDense->matrixL().transpose().solve(whitenoise);
      return 0;
    }
    
    outv.noalias() += _factorDense->matrixL() * whitenoise;
    return 0;
  }
  return 0;

}


MatrixSquareSymmetricSim::~MatrixSquareSymmetricSim()
{
  _clear();
}

void MatrixSquareSymmetricSim::_clear()
{
  delete _factorDense;
  delete _factorSparse;
  _factorDense = nullptr;
  _factorSparse = nullptr;
} 

void MatrixSquareSymmetricSim::_prepare() const //duplicated from MatrixSquareSymmetric
{
  if (isEmpty()) 
  {
    messerr("Your object is empty.");
    return;
  }
  if (!_sparse)
  {
    if (_factorDense == nullptr)
    {
      const Eigen::MatrixXd* a = ((AMatrixDense*)_mat)->getTab();
      _factorDense = new Eigen::LLT<Eigen::MatrixXd>();
      *_factorDense = a->llt();
    }
  }
  if(_sparse)
  {
    if (_factorSparse == nullptr)
    {
      _factorSparse = new Cholesky((MatrixSparse*)_mat);
    }
  }
}
