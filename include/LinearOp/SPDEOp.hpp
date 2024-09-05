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
  SPDEOp(const PrecisionOpMulti* const pop      = nullptr, 
         const ProjMulti*        const proj     = nullptr,
         const ALinearOp*        const invNoise = nullptr);
  virtual ~SPDEOp();

  int getSize() const override;
  VectorDouble kriging(const VectorDouble& dat) const;
  VectorDouble simulateCond(const VectorDouble& dat) const;

#ifndef SWIG
public:
  int kriging(const Eigen::VectorXd& inv,
                    Eigen::VectorXd& out) const;
protected:
  int _addToDest(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
private: 
  virtual int _solve(const Eigen::VectorXd& in,Eigen::VectorXd& out) const;
  int _buildRhs(const Eigen::VectorXd& inv) const;
#endif

private:
  void _prepare(bool w1 = true, bool w2 = true) const;
#ifndef SWIG
protected:
  virtual int _addToDestImpl(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const;
#endif

protected:
  const PrecisionOpMulti* const _Q;
  const ProjMulti*        const _Proj;
  const ALinearOp*        const _invNoise;

private:
  mutable Eigen::VectorXd _workdat1; 
  mutable Eigen::VectorXd _workdat2;
  mutable Eigen::VectorXd _rhs;

};

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(SPDEOp)
#endif