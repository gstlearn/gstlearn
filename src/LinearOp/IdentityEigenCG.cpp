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
#include <LinearOp/IdentityEigenCG.hpp>

IdentityEigenCG::IdentityEigenCG(int n) 
  : ALinearOpEigenCG<IdentityEigenCG>()
  , _n(n)
{
}

IdentityEigenCG::~IdentityEigenCG() 
{
}

/*****************************************************************************/
/*!
**  Evaluate the product (by the IdentityEigenCG) : 'outv' = I * 'inv' = 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
void IdentityEigenCG::_evalDirect(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const
{
  for(int i=0, n=_n; i<n; i++)
    outv[i] = 2*inv[i];
}

/* Force execution of Eigen conjugate gradient
void IdentityEigenCG::evalInverse(const VectorDouble &inv, VectorDouble &outv) const
{
  evalDirect(inv,outv);
}
*/