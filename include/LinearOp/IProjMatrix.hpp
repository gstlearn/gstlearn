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
#pragma once

#include "LinearOp/ALinearOp.hpp"
#include "gstlearn_export.hpp"
#include <span>

#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
  #include <Eigen/src/Core/Matrix.h>
  #include <Basic/VectorNumT.hpp>
#endif

class VectorEigen;
class GSTLEARN_EXPORT IProjMatrix
{
public:
  IProjMatrix() { }
  virtual ~IProjMatrix() { }
  virtual int point2mesh(const VectorDouble& inv, VectorDouble& outv) const; //TODO supress virtual
  virtual int mesh2point(const VectorDouble& inv, VectorDouble& outv) const; // when job is finished
  int point2mesh(const VectorEigen& inv, VectorEigen  & outv) const;
  int mesh2point(const VectorEigen& inv, VectorEigen   & outv)const; 
  int point2mesh(const constvect& inv,vect& out) const;
  int mesh2point(const constvect& inv,vect& out) const;

  virtual int getApexNumber() const = 0;
  virtual int getPointNumber() const = 0;

  #ifndef SWIG
  int mesh2point(const Eigen::VectorXd& inv,
                       Eigen::VectorXd& outv) const;
  int point2mesh(const Eigen::VectorXd& inv,
                       Eigen::VectorXd& outv) const;    
  int addMesh2point(const constvect& inv,
                    vect& outv) const;
  int addPoint2mesh(const constvect& inv,
                    vect& outv) const;       
 
  protected:
  virtual int _addPoint2mesh(const constvect& inv,
                             vect& outv) const
  {
    DECLARE_UNUSED(inv);
    DECLARE_UNUSED(outv);
    return 1;
  }
  virtual int _addMesh2point(const constvect& inv,
                             vect& outv) const
  {
    DECLARE_UNUSED(inv);
    DECLARE_UNUSED(outv);
    return 1;
  }
  #endif
};
