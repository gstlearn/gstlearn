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
#include <Eigen/src/Core/Map.h>

#ifndef SWIG
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
  LinearOpCGSolver(const TLinOP* linop);

  void solve(const VectorDouble& rhs, VectorDouble& out);
  void solve(const VectorEigen& rhs, VectorEigen& out);
  void setMaxIterations(int n) {cg.setMaxIterations(n);}
  void setTolerance(double tol) {cg.setTolerance(tol);}
  int  getIterations() const { return cg.iterations();}
  double getError() const { return  cg.error();}
#ifndef SWIG
  void solve(const Eigen::Map<const Eigen::VectorXd>& rhs,
             Eigen::Map<Eigen::VectorXd>& out);
  void solveWithGuess(const Eigen::Map<const Eigen::VectorXd>& rhs,
                      const Eigen::Map<const Eigen::VectorXd>& guess,
                      Eigen::Map<Eigen::VectorXd>& out);
private:
  Eigen::ConjugateGradient<TLinOP,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::IdentityPreconditioner> cg;
#endif
};

#ifndef SWIG
template<typename TLinOP>
LinearOpCGSolver<TLinOP>::LinearOpCGSolver(const TLinOP* linop)
{
  if (linop == nullptr)
    throw("linop must be valid and inherit from ALinearOpEigenCG to use Eigen CG");

  cg.compute(*linop);
}

template<typename TLinOP>
void LinearOpCGSolver<TLinOP>::solve(const VectorDouble& rhs, VectorDouble& out)
{
  Eigen::Map<const Eigen::VectorXd> myRhs(rhs.data(), rhs.size());
  Eigen::Map<Eigen::VectorXd> myOut(out.data(), out.size());
  // Assume outv has the good size
  solve(myRhs, myOut);
}

template<typename TLinOP>
void LinearOpCGSolver<TLinOP>::solve(const VectorEigen& rhs, VectorEigen& out)
{
  auto outmap = out.getMap();
  solve(rhs.getMap(), outmap);
}

template<typename TLinOP>
void LinearOpCGSolver<TLinOP>::solve(
  const Eigen::Map<const Eigen::VectorXd>& rhs,
  Eigen::Map<Eigen::VectorXd>& out)
{
  out = cg.solve(rhs);
}

template<typename TLinOP>
void LinearOpCGSolver<TLinOP>::solveWithGuess(const Eigen::Map<const Eigen::VectorXd>& rhs,
                                              const Eigen::Map<const Eigen::VectorXd>& guess,
                                              Eigen::Map<Eigen::VectorXd>& out)
{
  out = cg.solveWithGuess(rhs,guess);
}


#endif
