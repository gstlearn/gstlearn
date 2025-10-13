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

#ifndef SWIG
#include "LinearOp/ALinearOpEigenCG.hpp"
DECLARE_EIGEN_TRAITS(ScaleOp)
#else
#include "LinearOp/ALinearOp.hpp"
#endif

namespace gstlrn {

class GSTLEARN_EXPORT ScaleOp:
#ifndef SWIG
  public ALinearOpEigenCG<ScaleOp>
#else
  public ALinearOp
#endif
{

public:
  ScaleOp(Id n, double scale = 1.);
  virtual ~ScaleOp();

  Id getSize() const override { return _n; }

#ifndef SWIG
protected:
  Id _addToDest(const gstlrn::constvect inv, gstlrn::vect outv) const override;
#endif

private:
  Id _n;
  double _scale;
};

} // namespace gstlrn

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(ScaleOp)
#endif
