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

#include "LinearOp/ALinearOpEigenCG.hpp"

DECLARE_EIGEN_TRAITS(IdentityEigenCG)

class GSTLEARN_EXPORT IdentityEigenCG: public ALinearOpEigenCG<IdentityEigenCG>
{

public:
  IdentityEigenCG(int n);
  virtual ~IdentityEigenCG();

  //void evalInverse(const VectorDouble& inv, VectorDouble& outv) const override;
  int getSize() const override { return _n; }

protected:
  void _evalDirect(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;

private:
  int _n;
};

DECLARE_EIGEN_PRODUCT(IdentityEigenCG)