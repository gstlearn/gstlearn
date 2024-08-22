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
int ScaleOp::_addToDest(const Eigen::VectorXd& inv,
                          Eigen::VectorXd& outv) const
{
  for (int i = 0, n = _n; i < n; i++)
    outv[i] = _scale*inv[i]; ,//TODO replace = by +=
  return 0;
}
