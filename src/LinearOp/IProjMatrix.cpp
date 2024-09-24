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
#include "geoslib_define.h"

#include <Eigen/Core>
#include <Eigen/Dense>

int IProjMatrix::mesh2point(const VectorDouble& inv,
                                  VectorDouble& outv) const
{
  outv.resize(getPointNumber());
  constvect myInv(inv.data(), inv.size());
  vect myOut(outv);
  return mesh2point(myInv, myOut);  
}

int IProjMatrix::point2mesh(const VectorDouble& inv,
                           VectorDouble& outv) const
{
  outv.resize(getApexNumber());
  constvect myInv(inv.data(), inv.size());
  vect myOut(outv);
  return point2mesh(myInv, myOut); 
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
  outv.resize(getApexNumber());
  constvect invs(inv.data(),inv.size());
  vect outs(outv.data(),outv.size());
  return point2mesh(invs,outs);

}

int IProjMatrix::mesh2point(const Eigen::VectorXd& inv,
                                  Eigen::VectorXd& outv) const
{
  outv.resize(getPointNumber());
  constvect invs(inv.data(),inv.size());
  vect outvs(outv.data(),outv.size());
  return mesh2point(invs,outvs);
}

int IProjMatrix::addMesh2point(const constvect& inv,
                    vect& outv) const
{
  return _addMesh2point(inv,outv);
}

int IProjMatrix::addPoint2mesh(const constvect& inv,
                    vect& outv) const
{
  return _addPoint2mesh(inv, outv);
}    

int IProjMatrix::mesh2point(const constvect& inv,
                    vect& outv) const
{
  std::fill(outv.begin(),outv.end(),0.);
  return _addMesh2point(inv,outv);
}

int IProjMatrix::point2mesh(const constvect& inv,
                    vect& outv) const
{ 
  std::fill(outv.begin(),outv.end(),0.);
  return _addPoint2mesh(inv, outv);
}    