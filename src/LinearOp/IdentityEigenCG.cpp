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
#include "LinearOp/IdentityEigenCG.hpp"
#include "LinearOp/ALinearOpEigenCG.hpp"

IMPLEMENT_LINEAROP_EIGEN_CG_INTERFACE(IdentityEigenCG)

IdentityEigenCG::IdentityEigenCG(int n) 
  : ALinearOpEigenCG<IdentityEigenCG>()
  , _n(n)
{
}

IdentityEigenCG::~IdentityEigenCG() {}

/*****************************************************************************/
/*!
**  Evaluate the product (by the IdentityEigenCG) : 'outv' = I * 'inv' = 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
void IdentityEigenCG::_evalDirectEigen(const Eigen::VectorXd& inv,
                                       Eigen::VectorXd& outv) const
{
  /// TODO : Add a scale parameter and rename in ScaleEigenCG ?
  for (int i = 0, n = _n; i < n; i++)
    outv[i] = 2*inv[i];
}
