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
DECLARE_EIGEN_TRAITS(SPDEOp)
#else
#include "LinearOp/ALinearOp.hpp"
#endif

class GSTLEARN_EXPORT SPDEOp:
#ifndef SWIG
  public ALinearOpEigenCG<SPDEOp>
#else
  public ALinearOp
#endif
{

public:
  SPDEOp(int n, double scale = 1.);
  virtual ~SPDEOp();

  int getSize() const override { return _n; }

#ifndef SWIG
protected:
  int _addToDest(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
#endif

private:
  int _n;
  double _scale;
};

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(SPDEOp)
#endif