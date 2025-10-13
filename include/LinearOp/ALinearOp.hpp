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

#include "Basic/VectorNumT.hpp"
#include "geoslib_define.h"

#ifndef SWIG
#  include <Eigen/Core>
#  include <Eigen/src/Core/Matrix.h>
#endif

namespace gstlrn
{
class GSTLEARN_EXPORT ALinearOp
{
public:
  ALinearOp();
  ALinearOp(const ALinearOp& m)            = delete;
  ALinearOp& operator=(const ALinearOp& m) = delete;
  virtual ~ALinearOp() {}
  virtual Id getSize() const = 0;

  Id evalDirect(const VectorDouble& inv, VectorDouble& outv) const;
  VectorDouble evalDirect(const VectorDouble& in) const;
  virtual void multiplyByValueAndAddDiagonal(double v1 = 1., double v2 = 0.) const;
  virtual void resetModif() const;
  void setUseFactor(bool usefactor) { _usefactor = usefactor; }
#ifndef SWIG

public:
  Id evalDirect(constvect inv, vect outv) const;
  Id addToDest(const constvect inv, vect outv) const;
  Id addToDest(const ::Eigen::VectorXd& inv, ::Eigen::VectorXd& outv) const;

protected:
  virtual Id _addToDest(constvect inv, vect outv) const = 0;
#endif

private:
  mutable bool _usefactor;
  mutable double _idfactor;
  mutable double _factor;
  mutable VectorDouble _temp;
};
} // namespace gstlrn