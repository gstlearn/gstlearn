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

#include "Basic/VectorNumT.hpp"
#include "Matrix/VectorEigen.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <span>
#ifndef SWIG
#  include <Eigen/Core>
#  include <Eigen/Dense>
#endif

typedef const std::span<const double> constvect;
typedef std::span<double> vect ;
class GSTLEARN_EXPORT ALinearOp
{
public:
  virtual ~ALinearOp() {}
  virtual int getSize() const = 0;
  int evalDirect(constvect& inv, vect& outv) const;
  int evalDirect(const VectorDouble& inv, VectorDouble& outv) const;
  VectorDouble evalDirect(const VectorDouble& in) const;
  int evalDirect(const VectorEigen& inv, VectorEigen& outv) const;
  int addToDest(const VectorDouble& inv, VectorDouble& outv) const;
  int addToDest(const VectorEigen& inv, VectorEigen& outv) const;
#ifndef SWIG
  public:
  int evalDirect(const Eigen::VectorXd& inv,
                 Eigen::VectorXd& outv) const;
  int addToDest(const Eigen::VectorXd& inv,
                Eigen::VectorXd& outv) const;

protected:
  virtual int _addToDest(constvect& inv,
                         vect& outv) const = 0;
  virtual int _addToDest(const Eigen::VectorXd& inv,
                         Eigen::VectorXd & outv) const = 0;
#endif
};
