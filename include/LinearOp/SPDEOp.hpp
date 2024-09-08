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
#include "LinearOp/ASimulable.hpp"
#include "LinearOp/ISimulable.hpp"

#ifndef SWIG
#include "LinearOp/ALinearOpEigenCG.hpp"
DECLARE_EIGEN_TRAITS(SPDEOp)
#else
#include "LinearOp/ALinearOp.hpp"
#endif

#include "LinearOp/LinearOpCGSolver.hpp"


class PrecisionOpMulti;
class ProjMulti;


class GSTLEARN_EXPORT SPDEOp:
#ifndef SWIG
  public ALinearOpEigenCG<SPDEOp>, ASimulable
#else
  public ISimulable
#endif
{

public:
  SPDEOp(const PrecisionOpMulti* const pop      = nullptr, 
         const ProjMulti*        const proj     = nullptr,
         const ISimulable*       const invNoise = nullptr,
         bool  noiseToDelete                    = false);
  virtual ~SPDEOp();

  int getSize() const override;
  VectorDouble kriging(const VectorDouble& dat) const;
  VectorDouble krigingWithGuess(const VectorDouble& dat, const VectorDouble& guess) const;
  void setMaxIterations(int n) {_solver.setMaxIterations(n);}
  void setTolerance(double tol) {_solver.setTolerance(tol);}
  int  getIterations() const { return _solver.getIterations();}
  double getError() const { return  _solver.getError();}
#ifndef SWIG
public:
  int kriging(const Eigen::VectorXd& inv,
                    Eigen::VectorXd& out) const;
  int krigingWithGuess(const Eigen::VectorXd& inv,
                       const Eigen::VectorXd& guess,
                             Eigen::VectorXd& out) const;
protected:
  int _addToDest(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
  int _addSimulateToDest(const Eigen::VectorXd& whitenoise, Eigen::VectorXd& outv) const override;

private: 
  virtual int _solve(const Eigen::VectorXd& in,Eigen::VectorXd& out) const;
  int _solveWithGuess(const Eigen::VectorXd& in,const Eigen::VectorXd &guess,Eigen::VectorXd& out) const;

  int _buildRhs(const Eigen::VectorXd& inv) const;
#endif

private:
  void _prepare(bool w1 = true, bool w2 = true) const;
  void _fake() const override{};
#ifndef SWIG
private:
  virtual int _addToDestImpl(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const;
#endif

protected:
  const PrecisionOpMulti* const _Q;
  const ProjMulti*        const _Proj;
  const ISimulable*       const _invNoise;

private:
  bool    _noiseToDelete;
  mutable Eigen::VectorXd _workdat1; 
  mutable Eigen::VectorXd _workdat2;
  mutable Eigen::VectorXd _rhs;
  mutable LinearOpCGSolver<SPDEOp> _solver;

};

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(SPDEOp)
#endif