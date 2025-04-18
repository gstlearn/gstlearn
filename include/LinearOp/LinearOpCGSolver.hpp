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

// iostream is included here as it is used in Eigen function (std::cerr)
#include <iostream>

#ifndef SWIG
#  include <Eigen/Core>
#  include <Eigen/Dense>
#  include <Eigen/IterativeLinearSolvers>
#  include <Eigen/src/Core/Matrix.h>
#  include <unsupported/Eigen/IterativeSolvers>
#endif

class ALinearOpCGSolver
{
public:
  virtual ~ALinearOpCGSolver() = default;
  virtual void solve(const VectorDouble& rhs, VectorDouble& out) = 0;
  virtual void setMaxIterations(int n) = 0;
  virtual void setTolerance(double tol) = 0;
  virtual int  getIterations() const = 0;
  virtual double getError() const = 0;
#ifndef SWIG
  virtual void solve(const constvect in, const vect out) = 0;
  virtual void solve(const Eigen::Map<const Eigen::VectorXd>& rhs,
                     Eigen::Map<Eigen::VectorXd>& out) = 0;

  virtual void solveWithGuess(const constvect rhs, const constvect guess, vect out) = 0;
  virtual void solveWithGuess(const Eigen::Map<const Eigen::VectorXd>& rhs,
                              const Eigen::Map<const Eigen::VectorXd>& guess,
                              Eigen::Map<Eigen::VectorXd>& out) = 0;
#endif
};

template<typename TLinOP>
class LinearOpCGSolver : public ALinearOpCGSolver
{
public:
  LinearOpCGSolver(const TLinOP* linop);
  virtual ~LinearOpCGSolver() = default;

  void solve(const VectorDouble& rhs, VectorDouble& out) override;
  void setMaxIterations(int n) override {cg.setMaxIterations(n);}
  void setTolerance(double tol) override {cg.setTolerance(tol);}
  int  getIterations() const override { return cg.iterations();}
  double getError() const override { return  cg.error();}
#ifndef SWIG
  void solve(const constvect in, const vect out) override;
  void solve(const Eigen::Map<const Eigen::VectorXd>& rhs,
             Eigen::Map<Eigen::VectorXd>& out) override;

  void solveWithGuess(const constvect rhs, const constvect guess, vect out) override;
  void solveWithGuess(const Eigen::Map<const Eigen::VectorXd>& rhs,
                      const Eigen::Map<const Eigen::VectorXd>& guess,
                      Eigen::Map<Eigen::VectorXd>& out) override;
private:
  Eigen::ConjugateGradient<TLinOP,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::IdentityPreconditioner> cg;
#endif
};

#ifndef SWIG
template<typename TLinOP>
LinearOpCGSolver<TLinOP>::LinearOpCGSolver(const TLinOP* linop) : ALinearOpCGSolver()
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

template<typename TLinOP>
void LinearOpCGSolver<TLinOP>::solve(const constvect in, const vect out)
{
  Eigen::Map<const Eigen::VectorXd> inm(in.data(),in.size());
  Eigen::Map<Eigen::VectorXd> outm(out.data(),out.size());
  solve(inm,outm);
}

template<typename TLinOP>
void LinearOpCGSolver<TLinOP>::solveWithGuess(const constvect rhs,
                                              const constvect guess,
                                              vect out)
{
  Eigen::Map<const Eigen::VectorXd> rhsm(rhs.data(),rhs.size());
  Eigen::Map<const Eigen::VectorXd> guessm(guess.data(),guess.size());
  Eigen::Map<Eigen::VectorXd> outm(out.data(),out.size());
  outm = cg.solveWithGuess(rhsm,guessm);
}

#endif
