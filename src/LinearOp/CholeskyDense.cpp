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
  _prepare();
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
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(),vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor->matrixL().transpose().solve(mvecin);
  return 0;             
}

int CholeskyDense::addLtX(const constvect vecin, vect vecout) const
{
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor->matrixL().transpose() * mvecin;
  return 0;
}

int CholeskyDense::addLX(const constvect vecin, vect vecout) const
{
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(),vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor->matrixL() * mvecin;
  return 0;
}

int CholeskyDense::addInvLX(const constvect vecin, vect vecout) const
{
  Eigen::Map<const Eigen::VectorXd> mvecin(vecin.data(), vecin.size());
  Eigen::Map<Eigen::VectorXd> mvecout(vecout.data(), vecout.size());
  mvecout.noalias() += _factor->matrixL().solve(mvecin);
  return 0;
}

int CholeskyDense::addSolveX(const constvect vecin, vect vecout) const
{
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
  auto diag  = _factor->matrixLLT().diagonal();
  double det = 0.;
  for (int i = 0; i < _factor->rows(); i++) det += log(diag[i]);
  return 2. * det;
}

VectorDouble CholeskyDense::getCholeskyTL() const
{
  _computeTL();
  return _tl;
}

double CholeskyDense::getCholeskyTL(int i, int j) const
{
  _computeTL();
  int neq = _size;
  return (i >= j) ? _TL(i, j) : 0.;
}

double CholeskyDense::getCholeskyTL(int iad) const
{
  _computeTL();
  return _tl[iad];
}

VectorDouble CholeskyDense::getCholeskyXL() const
{
  _computeXL();
  return _xl;
}

double CholeskyDense::getCholeskyXL(int i, int j) const
{
  _computeXL();
  int neq = _size;
  return (i >= j) ? _XL(i, j) : 0.;
}

void CholeskyDense::_prepare() const
{
  const Eigen::MatrixXd* a = ((AMatrixDense*)_mat)->getTab();
  _factor                  = new Eigen::LLT<Eigen::MatrixXd>();
  *_factor                 = a->llt();
}

void CholeskyDense::_computeTL() const
{
  if (!_tl.empty()) return;
  int neq = _size;

  _tl.resize(_getTriangleSize());
  Eigen::MatrixXd mymat = _factor->matrixL();
  for (int ip = 0; ip < neq; ip++)
    for (int jp = 0; jp <= ip; jp++) _TL(ip, jp) = mymat(ip, jp);
}

void CholeskyDense::_computeXL() const
{
  if (!_xl.empty()) return;

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
 ** \param[in]  neq  number of equations in the system
 ** \param[in]  nrhs number of columns in x
 ** \param[in]  a    matrix (dimension neq * nrhs)
 **
 *****************************************************************************/
MatrixRectangular CholeskyDense::productCholeskyInPlace(int mode,
                                                        int neq,
                                                        int nrhs,
                                                        const MatrixRectangular& a)
{
  _computeTL();
  MatrixRectangular x(neq, nrhs);

  double val = 0.;
  if (mode == 0)
  {
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = i; j < neq; j++) val += _TL(j, i) * a.getValue(j, irhs);
        x.setValue(i, irhs, val);
      }
  }
  else if (mode == 1)
  {
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = 0; j <= i; j++) val += _TL(i, j) * a.getValue(j, irhs);
        x.setValue(i, irhs, val);
      }
  }
  else if (mode == 2)
  {
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = 0; j <= i; j++) val += a.getValue(irhs, j) * _TL(i, j);
        x.setValue(irhs, i, val);
      }
  }
  else if (mode == 3)
  {
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = i; j < neq; j++) val += a.getValue(irhs, j) * _TL(j, i);
        x.setValue(irhs, i, val);
      }
  }
  else if (mode == 4)
  {
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = 0; j <= i; j++) val += a.getValue(irhs, j) * _TL(i, j);
        x.setValue(irhs, i, val);
      }
  }
  else if (mode == 5)
  {
    for (int irhs = 0; irhs < nrhs; irhs++)
      for (int i = 0; i < neq; i++)
      {
        val = 0.;
        for (int j = i; j < neq; j++) val += a.getValue(irhs, j) * _TL(j, i);
        x.setValue(irhs, i, val);
      }
  }
  return x;
}

/*****************************************************************************/
/*!
 **  Performs the product B = TL * A * TU or TU * A * TL
 **  where TL,TU is a triangular matrix and A a square symmetric matrix
 **
 ** \param[in]  mode  0: TL * A * TU; 1: TU * A * TL
 ** \param[in]  neq  number of equations in the system
 ** \param[in]  a    Square symmetric matrix (optional)
 **
 *****************************************************************************/
MatrixSquareSymmetric CholeskyDense::normCholeskyInPlace(int mode,
                                                         int neq,
                                                         const MatrixSquareSymmetric& a)
{
  MatrixSquareSymmetric b(neq);
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
  return b;
}
