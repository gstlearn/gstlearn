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

ACholesky::ACholesky(const AMatrix* mat)
  : _mat(mat)
  , _size(0)
{
  _size = mat->getNRows();
}

int ACholesky::_addToDest(const constvect vecin, vect vecout) const
{
  return addLX(vecin, vecout);
}

int ACholesky::_addSimulateToDest(const constvect whitenoise, vect vecout) const
{
  return addInvLtX(whitenoise, vecout);
}

int ACholesky::solve(const constvect vecin, vect vecout) const
{
  std::fill(vecout.begin(), vecout.end(), 0.);
  return addSolveX(vecin, vecout);
}

int ACholesky::LX(const constvect whitenoise, vect vecout) const
{
  std::fill(vecout.begin(), vecout.end(), 0.);
  return addLX(whitenoise, vecout);
}

int ACholesky::InvLX(const constvect whitenoise, vect vecout) const
{
  std::fill(vecout.begin(), vecout.end(), 0.);
  return addInvLX(whitenoise, vecout);
}

int ACholesky::InvLtX(const constvect whitenoise, vect vecout) const
{
  std::fill(vecout.begin(), vecout.end(), 0.);
  return addInvLtX(whitenoise, vecout);
}

int ACholesky::LtX(const constvect whitenoise, vect vecout) const
{
  std::fill(vecout.begin(), vecout.end(), 0.);
  return addLtX(whitenoise, vecout);
}
