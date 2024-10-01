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

ACholesky::ACholesky(const AMatrix* mat, bool inverse)
  : _mat(mat)
  , _inverse(inverse)
  , _size(0)
{
  _size = mat->getNRows();
}

int ACholesky::_addToDest(const constvect vecin, vect vecout) const
{
  return _addLX(vecin, vecout);
}

int ACholesky::_addSimulateToDest(const constvect whitenoise, vect vecout) const
{
  if (_inverse) return _addInvLtX(whitenoise, vecout);
  return _addLX(whitenoise, vecout);
}

int ACholesky::solve(const constvect vecin, vect vecout) const
{
  return _addSolveX(vecin, vecout);
}
