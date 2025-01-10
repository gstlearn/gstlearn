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

#ifndef SWIG
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#endif

class GSTLEARN_EXPORT ALinearOp
{
public:
  ALinearOp();
  ALinearOp(const ALinearOp& op) = delete;
  virtual ~ALinearOp() {}
  virtual int getSize() const = 0;

  int evalDirect(const VectorDouble& inv, VectorDouble& outv) const;
  VectorDouble evalDirect(const VectorDouble& in) const;
  virtual void multiplyByValueAndAddDiagonal(double v1 = 1.,double v2 = 0.);
  virtual void resetModif();
  void setUseFactor(bool usefactor)
  {
    _usefactor = usefactor;
  }
#ifndef SWIG

public:
  int evalDirect(constvect inv, vect outv) const;
  int addToDest(const constvect inv, vect outv) const;
  int addToDest(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const;

protected:
  virtual int _addToDest(constvect inv, vect outv) const = 0;
#endif

private:
  mutable bool   _usefactor;
  mutable double _idfactor;
  mutable double _factor;
  mutable VectorDouble _temp;
};
