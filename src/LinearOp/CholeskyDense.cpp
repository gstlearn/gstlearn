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
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "geoslib_define.h"
#include <Eigen/src/Core/Matrix.h>

#define TRI(i)        (((i) * ((i) + 1)) / 2)
#define SQ(i, j, neq) ((j) * neq + (i))
#define _TL(i, j)     _tl[SQ(i, j, neq) - TRI(j)] /* for i >= j */
#define _XL(i, j)     _xl[SQ(i, j, neq) - TRI(j)] /* for i >= j */

namespace gstlrn
{
CholeskyDense::CholeskyDense(const MatrixSymmetric& mat)
  : ACholesky(mat)
  , _tl()
  , _xl()
  , _factor()
{
  (void)_prepare(mat);
}

CholeskyDense::CholeskyDense(const CholeskyDense& m)
  : ACholesky(m)
  , _tl(m._tl)
  , _xl(m._xl)
  , _factor(m._factor)
{
}

CholeskyDense& CholeskyDense::operator=(const CholeskyDense& m)
{
  if (this != &m)
  {
    ACholesky::operator=(m);
    _tl     = m._tl;
    _xl     = m._xl;
    _factor = m._factor;
    _empty  = m._empty;
  }
  return *this;
}

CholeskyDense::~CholeskyDense()
{
  _clear();
}

void CholeskyDense::clear()
{
  _clear();
}
void CholeskyDense::_clear()
{
  _empty = true;
}

bool CholeskyDense::empty() const
{
  return _empty;
}

MatrixDense CholeskyDense::inverse() const
{
  if (!isReady()) return MatrixDense();
  auto Id = Eigen::MatrixXd::Identity(_size, _size);
  MatrixDense res(_size, _size);
  res.eigenMat() = _factor.solve(Id);
  return res;
}
void CholeskyDense::solveMatInPlace(const MatrixDense& mat, MatrixDense& res) const
{
  if (!isReady()) return;
  res.eigenMat() = _factor.solve(mat.eigenMat());
}
Id CholeskyDense::addInvLtX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor.matrixL().transpose().solve(mvecin);
  return 0;
}

Id CholeskyDense::addLtX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor.matrixL().transpose() * mvecin;
  return 0;
}

Id CholeskyDense::addLX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor.matrixL() * mvecin;
  return 0;
}

Id CholeskyDense::addInvLX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor.matrixL().solve(mvecin);
  return 0;
}

Id CholeskyDense::addSolveX(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  Id size = static_cast<Id>(vecin.size());
  Eigen::Map<const Eigen::VectorXd> bm(vecin.data(), size);
  Eigen::Map<Eigen::VectorXd> xm(vecout.data(), size);
  xm += _factor.solve(bm);
  return 0;
}

Id CholeskyDense::_getTriangleSize() const
{
  Id neq = _size;
  Id tri = neq * (neq + 1) / 2;
  return tri;
}

double CholeskyDense::computeLogDeterminant() const
{
  if (!isReady()) return TEST;
  auto diag  = _factor.matrixLLT().diagonal();
  double det = 0.;
  for (Id i = 0; i < _factor.rows(); i++) det += log(diag[i]);
  return 2. * det;
}

VectorDouble CholeskyDense::getLowerTriangle() const
{
  if (_computeTL()) return VectorDouble();
  return _tl;
}

double CholeskyDense::getLowerTriangle(Id i, Id j) const
{
  if (_computeTL()) return TEST;
  Id neq = _size;
  return (i >= j) ? _TL(i, j) : 0.;
}

VectorDouble CholeskyDense::getUpperTriangleInverse() const
{
  if (_computeXL()) return VectorDouble();
  return _xl;
}

double CholeskyDense::getUpperTriangleInverse(Id i, Id j) const
{
  if (_computeXL()) return TEST;
  ;
  Id neq = _size;
  return (i >= j) ? _XL(i, j) : 0.;
}

Id CholeskyDense::_prepare(const MatrixSymmetric& mat) const
{
  const auto& a = mat.eigenMat();
  _factor       = a.llt();
  _setReady();
  return 0;
}

Id CholeskyDense::setMatrix(const MatrixSymmetric& mat)
{
  _size = mat.getNRows();
  return _prepare(mat);
}

Id CholeskyDense::_computeTL() const
{
  if (!_tl.empty()) return 0;
  if (!isReady()) return 1;
  Id neq = _size;

  _tl.resize(_getTriangleSize());
  Eigen::MatrixXd mymat = _factor.matrixL();
  for (Id ip = 0; ip < neq; ip++)
    for (Id jp = 0; jp <= ip; jp++) _TL(ip, jp) = mymat(ip, jp);
  return 0;
}

Id CholeskyDense::_computeXL() const
{
  if (!_xl.empty()) return 0;
  if (!isReady()) return 1;
  if (_computeTL()) return 1;

  Id neq = _size;
  _xl.resize(_getTriangleSize());

  for (Id i = 0; i < neq; i++)
  {
    for (Id j = 0; j < i; j++)
    {
      double sum = 0.;
      for (Id l = j; l < i; l++) sum += _TL(i, l) * _XL(l, j);
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
void CholeskyDense::matProductInPlace(Id mode,
                                      const MatrixDense& a,
                                      MatrixDense& x)
{
  if (_computeTL()) return;
  auto n1 = a.getNRows();
  auto n2 = a.getNCols();
  x.reset(n1, n2);

  Id neq;
  Id nrhs;
  double val = 0.;
  if (mode == 0)
  {
    neq  = n1;
    nrhs = n2;
    for (Id irhs = 0; irhs < nrhs; irhs++)
      for (Id i = 0; i < neq; i++)
      {
        val = 0.;
        for (Id j = i; j < neq; j++)
          val += _TL(j, i) * a.getValue(j, irhs);
        x.setValue(i, irhs, val);
      }
  }
  else if (mode == 1)
  {
    neq  = n1;
    nrhs = n2;
    for (Id irhs = 0; irhs < nrhs; irhs++)
      for (Id i = 0; i < neq; i++)
      {
        val = 0.;
        for (Id j = 0; j <= i; j++)
          val += _TL(i, j) * a.getValue(j, irhs);
        x.setValue(i, irhs, val);
      }
  }
  else if (mode == 2)
  {
    nrhs = n1;
    neq  = n2;
    for (Id irhs = 0; irhs < nrhs; irhs++)
      for (Id i = 0; i < neq; i++)
      {
        val = 0.;
        for (Id j = 0; j <= i; j++)
          val += a.getValue(irhs, j) * _TL(i, j);
        x.setValue(irhs, i, val);
      }
  }
  else if (mode == 3)
  {
    nrhs = n1;
    neq  = n2;
    for (Id irhs = 0; irhs < nrhs; irhs++)
      for (Id i = 0; i < neq; i++)
      {
        val = 0.;
        for (Id j = i; j < neq; j++)
          val += a.getValue(irhs, j) * _TL(j, i);
        x.setValue(irhs, i, val);
      }
  }
  else if (mode == 4)
  {
    nrhs = n1;
    neq  = n2;
    for (Id irhs = 0; irhs < nrhs; irhs++)
      for (Id i = 0; i < neq; i++)
      {
        val = 0.;
        for (Id j = 0; j <= i; j++)
          val += a.getValue(irhs, j) * _TL(i, j);
        x.setValue(irhs, i, val);
      }
  }
  else if (mode == 5)
  {
    nrhs = n1;
    neq  = n2;
    for (Id irhs = 0; irhs < nrhs; irhs++)
      for (Id i = 0; i < neq; i++)
      {
        val = 0.;
        for (Id j = i; j < neq; j++)
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
void CholeskyDense::normMatInPlace(Id mode,
                                   Id neq,
                                   const MatrixSymmetric& a,
                                   MatrixSymmetric& b)
{
  if (_computeTL()) return;
  b.resize(neq, neq);
  double vala;

  for (Id i = 0; i < neq; i++)
    for (Id j = 0; j < neq; j++)
    {
      double val = 0.;
      if (mode == 0)
      {
        for (Id l = 0; l <= j; l++)
          for (Id k = 0; k <= i; k++)
          {
            if (!a.empty())
              vala = a.getValue(k, l);
            else
              vala = static_cast<double>(k == l);
            val += _TL(i, k) * vala * _TL(j, l);
          }
      }
      else
      {
        for (Id l = j; l < neq; l++)
          for (Id k = i; k < neq; k++)
          {
            if (!a.empty())
              vala = a.getValue(k, l);
            else
              vala = static_cast<double>(k == l);
            val += _TL(k, i) * vala * _TL(l, j);
          }
      }
      b.setValue(i, j, val);
    }
}
} // namespace gstlrn