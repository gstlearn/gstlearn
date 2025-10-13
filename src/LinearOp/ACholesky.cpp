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
#include "LinearOp/ACholesky.hpp"
#include "Matrix/AMatrix.hpp"
#include "Matrix/MatrixDense.hpp"

namespace gstlrn{
ACholesky::ACholesky(const AMatrix& mat)
  : _size(0)
  , _ready(false)
{
  _size = mat.getNRows();
}

ACholesky::ACholesky(const ACholesky& m)
  : _size(m._size)
  , _ready(m._ready)
{
}

ACholesky& ACholesky::operator=(const ACholesky& m)
{
  if (this != &m)
  {
    _size  = m._size;
    _ready = m._ready;
  }
  return *this;
}

Id ACholesky::_addToDest(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  return addLX(vecin, vecout);
}

Id ACholesky::_addSimulateToDest(const constvect whitenoise, vect vecout) const
{
  if (!isReady()) return 1;
  return addInvLtX(whitenoise, vecout);
}

Id ACholesky::solve(const constvect vecin, vect vecout) const
{
  if (!isReady()) return 1;
  std::fill(vecout.begin(), vecout.end(), 0.);
  return addSolveX(vecin, vecout);
}

Id ACholesky::LX(const constvect whitenoise, vect vecout) const
{
  if (!isReady()) return 1;
  std::fill(vecout.begin(), vecout.end(), 0.);
  return addLX(whitenoise, vecout);
}

Id ACholesky::InvLX(const constvect whitenoise, vect vecout) const
{
  if (!isReady()) return 1;
  std::fill(vecout.begin(), vecout.end(), 0.);
  return addInvLX(whitenoise, vecout);
}

Id ACholesky::InvLtX(const constvect whitenoise, vect vecout) const
{
  if (!isReady()) return 1;
  std::fill(vecout.begin(), vecout.end(), 0.);
  return addInvLtX(whitenoise, vecout);
}

Id ACholesky::LtX(const constvect whitenoise, vect vecout) const
{
  if (!isReady()) return 1;
  std::fill(vecout.begin(), vecout.end(), 0.);
  return addLtX(whitenoise, vecout);
}

Id ACholesky::solveMatrix(const MatrixDense& b, MatrixDense& x) const
{
  if (!isReady()) return 1;

  auto nrows = b.getNRows();
  auto ncols = b.getNCols();
  x.resize(nrows, ncols);

  VectorDouble xcol(nrows);
  for (Id icol = 0; icol < ncols; icol++)
  {
    auto bcol = b.getViewOnColumn(icol);
    solve(bcol, xcol);
    x.setColumn(icol, xcol);
  }
  return 0;
}

VectorDouble ACholesky::invLtX(const VectorDouble& vecin) const
{
  constvect spin(vecin);
  VectorDouble vecout(_size, 0);
  vect spout(vecout);
  addInvLtX(spin, spout);
  return vecout;
}

VectorDouble ACholesky::LtX(const VectorDouble& vecin) const
{
  constvect spin(vecin);
  VectorDouble vecout(_size, 0);
  vect spout(vecout);
  addLtX(spin, spout);
  return vecout;
}

VectorDouble ACholesky::LX(const VectorDouble& vecin) const
{
  constvect spin(vecin);
  VectorDouble vecout(_size, 0);
  vect spout(vecout);
  addLX(spin, spout);
  return vecout;
}

VectorDouble ACholesky::invLX(const VectorDouble& vecin) const
{
  constvect spin(vecin);
  VectorDouble vecout(_size, 0);
  vect spout(vecout);
  addInvLX(spin, spout);
  return vecout;
}
VectorDouble ACholesky::solveX(const VectorDouble& vecin) const
{
  constvect spin(vecin);
  VectorDouble vecout(_size, 0);
  vect spout(vecout);
  addSolveX(spin, spout);
  return vecout;
}

double ACholesky::computeLogDet(Id nMC) const
{
  DECLARE_UNUSED(nMC);
  return computeLogDeterminant();
}
}