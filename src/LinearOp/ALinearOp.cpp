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

void ALinearOp::addToDest(const VectorDouble& inv, VectorDouble& outv) const
{
   try
  {
    Eigen::Map<const Eigen::VectorXd> myInv(inv.data(), inv.size());
    Eigen::VectorXd myOut(outv.size());
    
    // Assume outv has the good size
    _addToDest(myInv, myOut);
    
    Eigen::Map<Eigen::VectorXd>(outv.data(), outv.size()) = myOut;
  }
  catch (const std::string& str)
  {
    // TODO : Check if std::exception can be used
    messerr("%s", str.c_str());
  }
}

void ALinearOp::addToDest(const VectorEigen& inv, VectorEigen& outv) const
{
  _addToDest(inv.getVector(), outv.getVector());
}

void ALinearOp::evalDirect(const Eigen::VectorXd& inv,
                            Eigen::VectorXd& outv)
{
    for (int i=0;i<outv.size();i++)
    {
      outv[i] = 0.;
    }              
    _addToDest(inv,outv);      
}



void ALinearOp::addToDest(const Eigen::VectorXd& inv,
                        Eigen::VectorXd& outv) const
{
  _addToDest(inv,outv);
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
void ALinearOp::evalDirect(const VectorDouble& inv,
                           VectorDouble& outv)
{
  VectorHelper::fill(outv,0.,inv.size());
  addToDest(inv,outv);
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
void ALinearOp::evalDirect(const VectorEigen& inv,
                           VectorEigen& outv)
{
  evalDirect(inv.getVector(), outv.getVector());
}

