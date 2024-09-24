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
#include "LinearOp/ScaleOp.hpp"
#include "LinearOp/ALinearOp.hpp"

ScaleOp::ScaleOp(int n, double scale) :
  _n(n), _scale(scale)
{
}

ScaleOp::~ScaleOp() {}

/*****************************************************************************/
/*!
**  Evaluate the product (by the ScaleOp) : 'outv' += I * 'inv' = 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
int ScaleOp::_addToDest(const constvect& inv,
                        vect& outv) const
{
  for (int i = 0, n = _n; i < n; i++)
    outv[i] += _scale * inv[i];
  return 0;
}
