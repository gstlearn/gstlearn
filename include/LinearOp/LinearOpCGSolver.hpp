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

#include "Matrix/VectorEigen.hpp"
#include "gstlearn_export.hpp"

#include "LinearOp/ILinearOpEigenCG.hpp"

#ifndef SWIG
#  include "LinearOp/ALinearOpEigenCG.hpp"
#  include <Eigen/Core>
#  include <Eigen/Dense>
#  include <Eigen/IterativeLinearSolvers>
#  include <Eigen/src/Core/Matrix.h>
#  include <unsupported/Eigen/IterativeSolvers>
#endif

template<typename TLinOP>
class LinearOpCGSolver  // No Export because it's a template
{
public:
  LinearOpCGSolver(ILinearOpEigenCG* linop);

  void solve(const VectorDouble& rhs, VectorDouble& out);
  void solve(const VectorEigen& rhs, VectorEigen& out);

#ifndef SWIG
  void solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& out);

private:
  Eigen::ConjugateGradient<TLinOP,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::IdentityPreconditioner> cg;
#endif
};

#ifndef SWIG
template<typename TLinOP>
LinearOpCGSolver<TLinOP>::LinearOpCGSolver(ILinearOpEigenCG* linop)
{
  ALinearOpEigenCG<TLinOP>* op = dynamic_cast<ALinearOpEigenCG<TLinOP>*>(linop);
  if (op == nullptr)
    throw("linop must inherit from ALinearOpEigenCG to use Eigen CG");

  cg.compute(*op);
}

template<typename TLinOP>
void LinearOpCGSolver<TLinOP>::solve(const VectorDouble& rhs, VectorDouble& out)
{
  Eigen::Map<const Eigen::VectorXd> myRhs(rhs.data(), rhs.size());
  Eigen::VectorXd myOut;
  // Assume outv has the good size
  solve(myRhs, myOut);
  Eigen::Map<Eigen::VectorXd>(out.data(), out.size()) = myOut;
}

template<typename TLinOP>
void LinearOpCGSolver<TLinOP>::solve(const VectorEigen& rhs, VectorEigen& out)
{
  solve(rhs.getVector(), out.getVector());
}

template<typename TLinOP>
void LinearOpCGSolver<TLinOP>::solve(const Eigen::VectorXd& rhs,
                                     Eigen::VectorXd& out)
{
  out = cg.solve(rhs);
}

#endif
