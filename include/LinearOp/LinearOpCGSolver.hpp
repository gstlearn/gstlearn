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

#include "LinearOp/ALinearOp.hpp"

#ifndef SWIG
#  include "LinearOp/ALinearOpEigenCG.hpp"
#  include <Eigen/Core>
#  include <Eigen/Dense>
#  include <Eigen/IterativeLinearSolvers>
#  include <Eigen/src/Core/Matrix.h>
#  include <unsupported/Eigen/IterativeSolvers>
#endif

template<typename TLinOP>
class LinearOpCGSolver
{
public:
  LinearOpCGSolver(TLinOP* linop);

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
LinearOpCGSolver<TLinOP>::LinearOpCGSolver(TLinOP* linop)
{
  if (linop == nullptr)
    throw("linop must be valid and inherit from ALinearOpEigenCG to use Eigen CG");

  cg.compute(*linop);
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
