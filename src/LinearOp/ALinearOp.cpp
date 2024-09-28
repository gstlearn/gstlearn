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
#include "LinearOp/ALinearOp.hpp"
#include "Basic/VectorNumT.hpp"

VectorDouble  ALinearOp::evalDirect(const VectorDouble& in) const
{
  VectorDouble res;
  evalDirect(in,res);
  return res;
}

int ALinearOp::addToDest(const Eigen::VectorXd& inv,
                Eigen::VectorXd& outv) const
{
  constvect ins(inv.data(),inv.size());
  vect outs(outv.data(),outv.size());
  return addToDest(ins,outs);

}

int ALinearOp::addToDest(const constvect inv, vect outv) const
{
  return _addToDest(inv,outv);
}

int ALinearOp::evalDirect(constvect inv, vect outv) const
{
  std::fill(outv.begin(),outv.end(),0.);
  return addToDest(inv, outv);
}

/*****************************************************************************/
/*!
**  Evaluate the product: 'outv' = Q * 'inv'
**
** \param[in]  inv     Array of input values
**
** \param[out] outv    Array of output values
**
*****************************************************************************/
int ALinearOp::evalDirect(const VectorDouble& inv,
                           VectorDouble& outv) const
{ 
  outv.resize(inv.size());
  constvect in(inv);
  vect out(outv);
  return evalDirect(in,out);
}


