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

#include "LinearOp/IProjMatrix.hpp"

#include "Matrix/VectorEigen.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

int IProjMatrix::mesh2point(const VectorDouble& inv,
                           VectorDouble& outv) const
{
    Eigen::Map<const Eigen::VectorXd> myInv(inv.data(), inv.size());
    Eigen::VectorXd myOut(outv.size());
    
    // Assume outv has the good size
    int error = _mesh2point(myInv, myOut);
    
    Eigen::Map<Eigen::VectorXd>(outv.data(), outv.size()) = myOut;
    return error;
}

int IProjMatrix::point2mesh(const VectorDouble& inv,
                           VectorDouble& outv) const
{
    Eigen::Map<const Eigen::VectorXd> myInv(inv.data(), inv.size());
    Eigen::VectorXd myOut(outv.size());
    
    // Assume outv has the good size
    int error=_point2mesh(myInv, myOut);
    
    Eigen::Map<Eigen::VectorXd>(outv.data(), outv.size()) = myOut;
    return error;
}

int IProjMatrix::mesh2point(const VectorEigen& inv, VectorEigen& outv) const
{
    return mesh2point(inv.getVector(), outv.getVector());
}

int IProjMatrix::point2mesh(const VectorEigen& inv, VectorEigen& outv) const
{
    return point2mesh(inv.getVector(), outv.getVector());
}

int IProjMatrix::point2mesh(const Eigen::VectorXd& inv,
                                  Eigen::VectorXd& outv) const
{
  return _point2mesh(inv,outv);

}

int IProjMatrix::mesh2point(const Eigen::VectorXd& inv,
                                  Eigen::VectorXd& outv) const
{
  return _mesh2point(inv,outv);
}