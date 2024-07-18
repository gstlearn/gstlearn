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

IdentityEigenCG::IdentityEigenCG(int n) 
  : ALinearOpEigenCG<IdentityEigenCG>()
  , _n(n)
{
}

IdentityEigenCG::~IdentityEigenCG() 
{
}
void IdentityEigenCG::evalInverse(const VectorDouble& inv,
                                  VectorDouble& outv) const
{
  ALinearOpEigenCG<IdentityEigenCG>::evalInverse(inv, outv);
}

void IdentityEigenCG::evalInverse(const VectorEigen& inv,
                                  VectorEigen& outv) const
{
  ALinearOpEigenCG<IdentityEigenCG>::evalInverse(inv, outv);
}

void IdentityEigenCG::evalDirect(const VectorDouble& inv,
                                 VectorDouble& outv) const
{
  ALinearOpEigenCG<IdentityEigenCG>::evalDirect(inv, outv);
}

void IdentityEigenCG::evalDirect(const VectorEigen& inv,
                                 VectorEigen& outv) const
{
  ALinearOpEigenCG<IdentityEigenCG>::evalDirect(inv, outv);
}

void IdentityEigenCG::setX0(const VectorDouble& x0)
{
  ALinearOpEigenCG<IdentityEigenCG>::setX0(x0);
}

void IdentityEigenCG::mustShowStats(bool status)
{
  ALinearOpEigenCG<IdentityEigenCG>::mustShowStats(status);
}

const LogStats& IdentityEigenCG::getLogStats() const
{
  return ALinearOpEigenCG<IdentityEigenCG>::getLogStats();
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
void IdentityEigenCG::_evalDirectEigen(const Eigen::VectorXd& inv,
                                        Eigen::VectorXd& outv) const
{
for (int i = 0, n = _n; i < n; i++)
  /// TODO : Add a scale parameter and rename in ScaleEigenCG ?
  outv[i] = 2*inv[i];
}
