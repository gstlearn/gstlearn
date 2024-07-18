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
class GSTLEARN_EXPORT LinearOpEigenCGSolver
{
public:
  LinearOpEigenCGSolver(ILinearOpEigenCG* linop);

  void solve(const VectorEigen& rhs, VectorEigen& out);

private:
  Eigen::ConjugateGradient<TLinOP,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::IdentityPreconditioner> cg;
};

#ifndef SWIG
template<typename TLinOP>
LinearOpEigenCGSolver<TLinOP>::LinearOpEigenCGSolver(ILinearOpEigenCG* linop)
{
  ALinearOpEigenCG<TLinOP>* op = dynamic_cast<ALinearOpEigenCG<TLinOP>*>(linop);
  cg.compute(*op);
}

template<typename TLinOP>
void LinearOpEigenCGSolver<TLinOP>::solve(const VectorEigen& rhs, VectorEigen& out)
{
  out.getVector() = cg.solve(rhs.getVector());
}
#endif
