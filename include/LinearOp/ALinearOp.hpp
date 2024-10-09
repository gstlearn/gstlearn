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

#include "geoslib_define.h"
#include "Basic/VectorNumT.hpp"

#ifndef swig
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#endif

class GSTLEARN_EXPORT ALinearOp
{
public:
  virtual ~ALinearOp() {}
  virtual int getSize() const = 0;

  int evalDirect(const VectorDouble& inv, VectorDouble& outv) const;
  VectorDouble evalDirect(const VectorDouble& in) const;

#ifndef SWIG

public:
  int evalDirect(constvect inv, vect outv) const;
  int addToDest(const constvect inv, vect outv) const;
  int addToDest(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const;

protected:
  virtual int _addToDest(constvect inv, vect outv) const = 0;
#endif
};
