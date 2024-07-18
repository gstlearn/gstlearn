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
DECLARE_EIGEN_TRAITS(IdentityEigenCG)
#else
#include "LinearOp/ILinearOpEigenCG.hpp"
#endif

class GSTLEARN_EXPORT IdentityEigenCG:
#ifndef SWIG
  public ALinearOpEigenCG<IdentityEigenCG>
#else
  public ILinearOpEigenCG
#endif
{

public:
  IdentityEigenCG(int n);
  virtual ~IdentityEigenCG();

  int getSize() const override { return _n; }

  void evalInverse(const VectorDouble& inv, VectorDouble& outv) const override;
  void evalInverse(const VectorEigen& inv, VectorEigen& outv) const override;

  void evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;
  void evalDirect(const VectorEigen& inv, VectorEigen& outv) const override;

  void setX0(const VectorDouble& x0) override;
  void mustShowStats(bool status) override;

  const LogStats& getLogStats() const override;

#ifndef SWIG
protected:
  void _evalDirectEigen(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
#endif

private:
  int _n;
};

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(IdentityEigenCG)
#endif