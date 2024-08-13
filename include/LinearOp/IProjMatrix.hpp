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

#include "gstlearn_export.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <Basic/VectorNumT.hpp>

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
  virtual int getApexNumber() const = 0;
  virtual int getPointNumber() const = 0;

  #ifndef SWIG
  int mesh2point(const Eigen::VectorXd& inv,
                       Eigen::VectorXd& outv) const;
  int point2mesh(const Eigen::VectorXd& inv,
                       Eigen::VectorXd& outv) const;             
  protected:
  virtual int _point2mesh(const Eigen::VectorXd& inv,
                                Eigen::VectorXd& outv) const = 0;
  virtual int _mesh2point(const Eigen::VectorXd& inv,
                                Eigen::VectorXd& outv) const = 0;
  #endif
};
