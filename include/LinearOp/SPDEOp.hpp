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

#include "LinearOp/IProjMatrix.hpp"
#include "LinearOp/ProjMulti.hpp"
#include "gstlearn_export.hpp"

#ifndef SWIG
#include "LinearOp/ALinearOpEigenCG.hpp"
DECLARE_EIGEN_TRAITS(SPDEOp)
#else
#include "LinearOp/ALinearOp.hpp"
#endif

class PrecisionOpMulti;
class ProjMulti;


class GSTLEARN_EXPORT SPDEOp:
#ifndef SWIG
  public ALinearOpEigenCG<SPDEOp>
#else
  public ALinearOp
#endif
{

public:
  SPDEOp(const PrecisionOpMulti* pop = nullptr, const ProjMulti* A = nullptr, const ALinearOp* invNoise = nullptr);
  virtual ~SPDEOp();

  int getSize() const override;

#ifndef SWIG
protected:
  virtual int _addToDest(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
#endif

protected:
  const PrecisionOpMulti* _Q;
  const ProjMulti*        _A;
  const ALinearOp*        _invNoise;
};

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(SPDEOp)
#endif