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
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Matrix/VectorEigen.hpp"

VectorDouble  ALinearOp::evalDirect(const VectorDouble& in)
{
  VectorDouble res(in.size());
  evalDirect(in,res);
  return res;
}

int ALinearOp::addToDest(const VectorDouble& inv, VectorDouble& outv) const
{
   try
  {
    Eigen::Map<const Eigen::VectorXd> myInv(inv.data(), inv.size());
    Eigen::VectorXd myOut(outv.size());
    VectorEigen::fill(myOut, 0.);
    // Assume outv has the good size
    if(_addToDest(myInv, myOut))
      return 1;
    
    VectorEigen::copy(myOut,outv);
  }
  catch (const std::string& str)
  {
    // TODO : Check if std::exception can be used
    messerr("%s", str.c_str());
  }
  return 0;
}

int ALinearOp::addToDest(const VectorEigen& inv, VectorEigen& outv) const
{
  return _addToDest(inv.getVector(), outv.getVector());
}

int ALinearOp::evalDirect(const Eigen::VectorXd& inv,
                            Eigen::VectorXd& outv)
{
    for (int i=0;i<(int)outv.size();i++)
    {
      outv[i] = 0.;
    }              
    return _addToDest(inv,outv);      
}



int ALinearOp::addToDest(const Eigen::VectorXd& inv,
                        Eigen::VectorXd& outv) const
{
  return _addToDest(inv,outv);
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
                           VectorDouble& outv)
{
  return addToDest(inv,outv);
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
int ALinearOp::evalDirect(const VectorEigen& inv,
                           VectorEigen& outv)
{
  return evalDirect(inv.getVector(), outv.getVector());
}

