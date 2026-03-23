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

#include "Basic/Message.hpp"
#include "Basic/VectorNumT.hpp"
#include "LinearOp/IdentityPrecond.hpp"

// iostream is included here as it is used in Eigen function (std::cerr)
#include <iostream>

#ifndef SWIG
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h>
#include <unsupported/Eigen/IterativeSolvers>

#endif

namespace gstlrn
{

  class ALinearOpCGSolver
  {
  public:
    virtual ~ALinearOpCGSolver() = default;
    virtual void solve(const VectorDouble& rhs, VectorDouble& out) = 0;
    virtual void setMaxIterations(Id n) = 0;
    virtual void setTolerance(double tol) = 0;
    virtual Id getIterations() const = 0;
    virtual double getError() const = 0;

    void setVerbose(bool v) { _verbose = v; }
#ifndef SWIG
    virtual void solve(const constvect in, const vect out) = 0;
    virtual void solve(const Eigen::Map<const Eigen::VectorXd>& rhs,
                       Eigen::Map<Eigen::VectorXd>& out) = 0;

    virtual void
      solveWithGuess(const constvect rhs, const constvect guess, vect out) = 0;
    virtual void solveWithGuess(const Eigen::Map<const Eigen::VectorXd>& rhs,
                                const Eigen::Map<const Eigen::VectorXd>& guess,
                                Eigen::Map<Eigen::VectorXd>& out) = 0;
    bool _verbose;
#endif
  };

  template<typename TLinOP, typename Precond = IdentityPreconditioner>
  class LinearOpCGSolver: public ALinearOpCGSolver
  {
  public:
    LinearOpCGSolver(const TLinOP* linop);
    virtual ~LinearOpCGSolver() = default;

    void solve(const VectorDouble& rhs, VectorDouble& out) override;

    void setMaxIterations(Id n) override { _cg.setMaxIterations(n); }

    void setTolerance(double tol) override { _cg.setTolerance(tol); }

    Id getIterations() const override { return _cg.iterations(); }

    Precond& getPreconditioner() { return _cg.preconditioner(); }

    double getError() const override { return _cg.error(); }
#ifndef SWIG
    void solve(const constvect in, const vect out) override;
    void solve(const Eigen::Map<const Eigen::VectorXd>& rhs,
               Eigen::Map<Eigen::VectorXd>& out) override;

    void solveWithGuess(const constvect rhs,
                        const constvect guess,
                        vect out) override;
    void solveWithGuess(const Eigen::Map<const Eigen::VectorXd>& rhs,
                        const Eigen::Map<const Eigen::VectorXd>& guess,
                        Eigen::Map<Eigen::VectorXd>& out) override;

  private:
    Eigen::ConjugateGradient<TLinOP, Eigen::Lower | Eigen::Upper, Precond> _cg;
#endif
  };
} // namespace gstlrn

using namespace gstlrn;
#ifndef SWIG
template<typename TLinOP, typename Precond>
LinearOpCGSolver<TLinOP, Precond>::LinearOpCGSolver(const TLinOP* linop)
  : ALinearOpCGSolver()
{
  if (linop == nullptr)
    throw(
      "linop must be valid and inherit from ALinearOpEigenCG to use Eigen CG");

  _cg.compute(*linop);
}

template<typename TLinOP, typename Precond>
void LinearOpCGSolver<TLinOP, Precond>::solve(const VectorDouble& rhs,
                                              VectorDouble& out)
{
  ::Eigen::Map<const ::Eigen::VectorXd> myRhs(rhs.data(), rhs.size());
  ::Eigen::Map<::Eigen::VectorXd> myOut(out.data(), out.size());
  // Assume outv has the good size
  solve(myRhs, myOut);
}

template<typename TLinOP, typename Precond>
void LinearOpCGSolver<TLinOP, Precond>::solve(
  const ::Eigen::Map<const ::Eigen::VectorXd>& rhs,
  ::Eigen::Map<Eigen::VectorXd>& out)
{
  out = _cg.solve(rhs);
  if (_verbose)
  {
    mestitle(0, "LinearOpCGSolver");
    messageFlush("Number of iterations: " + std::to_string(_cg.iterations())
                 + "/" + std::to_string(_cg.maxIterations()) + "\n");
    messageFlush("Error: " + std::to_string(_cg.error()) + "\n");
  }
}

template<typename TLinOP, typename Precond>
void LinearOpCGSolver<TLinOP, Precond>::solveWithGuess(
  const Eigen::Map<const Eigen::VectorXd>& rhs,
  const Eigen::Map<const Eigen::VectorXd>& guess,
  Eigen::Map<Eigen::VectorXd>& out)
{
  out = _cg.solveWithGuess(rhs, guess);
}

template<typename TLinOP, typename Precond>
void
  LinearOpCGSolver<TLinOP, Precond>::solve(const constvect in, const vect out)
{
  Eigen::Map<const Eigen::VectorXd> inm(in.data(), in.size());
  Eigen::Map<Eigen::VectorXd> outm(out.data(), out.size());
  solve(inm, outm);
}

template<typename TLinOP, typename Precond>
void LinearOpCGSolver<TLinOP, Precond>::solveWithGuess(const constvect rhs,
                                                       const constvect guess,
                                                       vect out)
{
  ::Eigen::Map<const Eigen::VectorXd> rhsm(rhs.data(), rhs.size());
  ::Eigen::Map<const Eigen::VectorXd> guessm(guess.data(), guess.size());
  ::Eigen::Map<Eigen::VectorXd> outm(out.data(), out.size());
  outm = _cg.solveWithGuess(rhsm, guessm);
}

#endif
