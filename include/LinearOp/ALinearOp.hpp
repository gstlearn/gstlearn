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

#ifndef SWIG
#  include <Eigen/Core>
#  include <Eigen/Dense>
#endif

class GSTLEARN_EXPORT ALinearOp
{
public:
  virtual ~ALinearOp() {}
  virtual int getSize() const = 0;
  
  virtual void evalDirect(const VectorDouble& inv, VectorDouble& outv) const;
  virtual void evalDirect(const VectorEigen& inv, VectorEigen& outv) const;

#ifndef SWIG
  virtual void evalDirect(const Eigen::VectorXd& inv,
                          Eigen::VectorXd& outv) const;

protected:
  virtual void _evalDirect(const Eigen::VectorXd& inv,
                           Eigen::VectorXd& outv) const = 0;
#endif
};
