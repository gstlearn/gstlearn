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
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "geoslib_define.h"
#include <Eigen/src/Core/Matrix.h>

#define TRI(i)        (((i) * ((i) + 1)) / 2)
#define SQ(i, j, neq) ((j)*neq + (i))
#define _TL(i, j) _tl[SQ(i, j, neq) - TRI(j)] /* for i >= j */
#define _XL(i, j) _xl[SQ(i, j, neq) - TRI(j)] /* for i >= j */

CholeskyDense::CholeskyDense(const MatrixSquareSymmetric* mat)
  : ACholesky(mat)
  , _tl()
  , _xl()
  , _factor(nullptr)
{
  (void) _prepare();
}

CholeskyDense::CholeskyDense(const CholeskyDense& m)
  : ACholesky(m)
  , _tl(m._tl)
  , _xl(m._xl)
  , _factor()
{
  if (m._factor != nullptr)
  {
    _factor = new Eigen::LLT<Eigen::MatrixXd>();
    _factor = m._factor;
  }
}

CholeskyDense& CholeskyDense::operator=(const CholeskyDense& m)
{
  if (this != &m)
  {
    ACholesky::operator=(m);
    _tl   = m._tl;
    _xl   = m._xl;
    if (m._factor != nullptr)
    {
      _factor = new Eigen::LLT<Eigen::MatrixXd>();
      _factor = m._factor;
    }
  }
  return *this;
}

CholeskyDense::~CholeskyDense()
{
  _clear();
}

void CholeskyDense::_clear()
{
  delete _factor;
  _factor = nullptr;
}

int CholeskyDense::addInvLtX(const constvect vecin, vect vecout) const
{
  if (! isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(),vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor->matrixL().transpose().solve(mvecin);
  return 0;             
}

int CholeskyDense::addLtX(const constvect vecin, vect vecout) const
{
  if (! isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor->matrixL().transpose() * mvecin;
  return 0;
}

int CholeskyDense::addLX(const constvect vecin, vect vecout) const
{
  if (! isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(),vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor->matrixL() * mvecin;
  return 0;
}

int CholeskyDense::addInvLX(const constvect vecin, vect vecout) const
{
  if (! isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor->matrixL().solve(mvecin);
  return 0;
}

int CholeskyDense::addSolveX(const constvect vecin, vect vecout) const
{
  if (! isReady()) return 1;
  int size = (int)vecin.size();
  Eigen::Map<const Eigen::VectorXd> bm(vecin.data(), size);
  Eigen::Map<Eigen::VectorXd> xm(vecout.data(), size);
  xm += _factor->solve(bm);
  return 0;
}

int CholeskyDense::_getTriangleSize() const
{
  int neq  = _size;
  int tri = neq * (neq + 1) / 2;
  return tri;
}

double CholeskyDense::computeLogDeterminant() const
{
  if (! isReady()) return TEST;
  auto diag  = _factor->matrixLLT().diagonal();
  double det = 0.;
  for (int i = 0; i < _factor->rows(); i++) det += log(diag[i]);
  return 2. * det;
}

VectorDouble CholeskyDense::getLowerTriangle() const
{
  if (_computeTL()) return VectorDouble();
  return _tl;
}

double CholeskyDense::getLowerTriangle(int i, int j) const
{
  if (_computeTL()) return TEST;
  int neq = _size;
  return (i >= j) ? _TL(i, j) : 0.;
}

VectorDouble CholeskyDense::getUpperTriangleInverse() const
{
  if (_computeXL()) return VectorDouble();
  return _xl;
}

double CholeskyDense::getUpperTriangleInverse(int i, int j) const
{
  if (_computeXL()) return TEST;;
  int neq = _size;
  return (i >= j) ? _XL(i, j) : 0.;
}

int CholeskyDense::_prepare() const
{
  if (_mat == nullptr) return 1;
  const Eigen::MatrixXd* a = ((AMatrixDense*)_mat)->getTab();
  _factor                  = new Eigen::LLT<Eigen::MatrixXd>();
  *_factor                 = a->llt();
  if (_factor == nullptr)
  {
    messerr("Error when computing the Cholesky Decmposition");
    return 1;
  }
  _setReady();
  return 0;
}

int CholeskyDense::setMatrix(const MatrixSquareSymmetric* mat)
{
  _mat = mat;
  _size = mat->getNRows();
  return _prepare();
}

int CholeskyDense::_computeTL() const
{
  if (!_tl.empty()) return 0;
  if (! isReady()) return 1;
  int neq = _size;

  _tl.resize(_getTriangleSize());
  Eigen::MatrixXd mymat = _factor->matrixL();
  for (int ip = 0; ip < neq; ip++)
    for (int jp = 0; jp <= ip; jp++) _TL(ip, jp) = mymat(ip, jp);
  return 0;
}

int CholeskyDense::_computeXL() const
{
  if (!_xl.empty()) return 0;
  if (!isReady()) return 1;
  if (_computeTL()) return 1;
  
  int neq = _size;
  _xl.resize(_getTriangleSize());

  for (int i = 0; i < neq; i++)
  {
    for (int j = 0; j < i; j++)
    {
      double sum = 0.;
      for (int l = j; l < i; l++) sum += _TL(i, l) * _XL(l, j);
      _XL(i, j) = -sum / _TL(i, i);
    }
    _XL(i, i) = 1. / _TL(i, i);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Performs the product between a triangular and a square matrix
 **  TL is the lower triangular matrix and X is a square matrix
 **
 ** \param[in]  mode Type of calculations:
 **             0 : X=TU%*%A
 **             1 : X=TL%*%A
 **             2 : X=A%*%TU
 **             3 : X=A%*%TL
 **             4 : X=t(A)%*%TU
 **             5 : X=t(A)%*%TL
 ** \param[in]  a    Input matrix
 ** \param[out] x    Resulting matrix (resized if necessary)
 **
 ** \remark The dimensions of 'a' and 'x' must match
 ** \remark Anyhow 'x' is resized to the same dimension as 'a'
 ** 
 *****************************************************************************/
void CholeskyDense::matProductInPlace(int mode,
                                      const MatrixRectangular& a,
                                      MatrixRectangular& x)
{
  if (_computeTL()) return;
  int n1 = a.getNRows();
  int n2 = a.getNCols();
  x.reset(n1, n2);

  int neq;
  int nrhs;
  double val = 0.;
  if (mode == 0)
  {
    neq  = n1;
    nrhs = n2;
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = i; j < neq; j++)
          val += _TL(j, i) * a.getValue(j, irhs);
        x.setValue(i, irhs, val);
      }
  }
  else if (mode == 1)
  {
    neq  = n1;
    nrhs = n2;
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = 0; j <= i; j++)
          val += _TL(i, j) * a.getValue(j, irhs);
        x.setValue(i, irhs, val);
      }
  }
  else if (mode == 2)
  {
    nrhs = n1;
    neq  = n2;
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = 0; j <= i; j++)
          val += a.getValue(irhs, j) * _TL(i, j);
        x.setValue(irhs, i, val);
      }
  }
  else if (mode == 3)
  {
    nrhs = n1;
    neq  = n2;
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = i; j < neq; j++)
          val += a.getValue(irhs, j) * _TL(j, i);
        x.setValue(irhs, i, val);
      }
  }
  else if (mode == 4)
  {
    nrhs = n1;
    neq  = n2;
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = 0; j <= i; j++)
          val += a.getValue(irhs, j) * _TL(i, j);
        x.setValue(irhs, i, val);
      }
  }
  else if (mode == 5)
  {
    nrhs = n1;
    neq  = n2;
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = i; j < neq; j++)
          val += a.getValue(irhs, j) * _TL(j, i);
        x.setValue(irhs, i, val);
      }
  }
}

/*****************************************************************************/
/*!
 **  Performs the product B = TL * A * TU or TU * A * TL
 **  where TL,TU is a triangular matrix and A a square symmetric matrix
 **
 ** \param[in]  mode  0: TL * A * TU; 1: TU * A * TL
 ** \param[in]  neq  number of equations in the system
 ** \param[in]  a    Square symmetric matrix (optional)
 ** \param[out] b    Square output matrix (resized if needed)
 **
 ** \remark The value of 'neq' could be derived from the input matrix 'a'
 ** \remark but this matrix is optional, hence the presence of argument 'neq' 
 **
 *****************************************************************************/
void CholeskyDense::normMatInPlace(int mode,
                                        int neq,
                                        const MatrixSquareSymmetric& a,
                                        MatrixSquareSymmetric& b)
{
  if (_computeTL()) return;
  b.resize(neq, neq);
  double vala;

  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      double val = 0.;
      if (mode == 0)
      {
        for (int l = 0; l <= j; l++)
          for (int k = 0; k <= i; k++)
          {
            if (!a.empty())
              vala = a.getValue(k, l);
            else
              vala = (double)(k == l);
            val += _TL(i, k) * vala * _TL(j, l);
          }
      }
      else
      {
        for (int l = j; l < neq; l++)
          for (int k = i; k < neq; k++)
          {
            if (!a.empty())
              vala = a.getValue(k, l);
            else
              vala = (double)(k == l);
            val += _TL(k, i) * vala * _TL(l, j);
          }
      }
      b.setValue(i, j, val);
    }
}

