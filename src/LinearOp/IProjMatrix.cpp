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

#include "geoslib_define.h"

int IProjMatrix::mesh2point(const VectorDouble& inv,
                                  VectorDouble& outv) const
{
  outv.resize(getPointNumber());
  return mesh2point(inv.getVector(), outv.getVector());
}

int IProjMatrix::point2mesh(const VectorDouble& inv,
                           VectorDouble& outv) const
{
  outv.resize(getApexNumber());
  return point2mesh(inv.getVector(), outv.getVector());
}

int IProjMatrix::addMesh2point(const constvect inv, vect outv) const
{
  return _addMesh2point(inv,outv);
}

int IProjMatrix::addPoint2mesh(const constvect inv, vect outv) const
{
  return _addPoint2mesh(inv, outv);
}

int IProjMatrix::mesh2point(const constvect inv, vect outv) const
{
  std::fill(outv.begin(),outv.end(),0.);
  return _addMesh2point(inv,outv);
}

int IProjMatrix::point2mesh(const constvect inv, vect outv) const
{ 
  std::fill(outv.begin(),outv.end(),0.);
  return _addPoint2mesh(inv, outv);
}
