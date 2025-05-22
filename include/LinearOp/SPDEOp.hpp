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

#include "LinearOp/ALinearOpEigenCG.hpp"
#include "Matrix/MatrixDense.hpp"

#include "Basic/VectorNumT.hpp"
#include "LinearOp/ASimulable.hpp"

#ifndef SWIG
#  include "LinearOp/ASimulableEigenCG.hpp"
DECLARE_EIGEN_TRAITS(SPDEOp)
#endif

#include "LinearOp/LinearOpCGSolver.hpp"

class PrecisionOpMulti;
class ProjMulti;


class GSTLEARN_EXPORT ASPDEOp : public virtual ALinearOp
{

public:
  ASPDEOp(const PrecisionOpMulti* const popKriging = nullptr,
          const ProjMulti* const projKriging       = nullptr,
          const ASimulable* const invNoise         = nullptr,
          const PrecisionOpMulti* const popSimu    = nullptr,
          const ProjMulti* const projSimu          = nullptr,
          bool noiseToDelete                       = false);
  virtual ~ASPDEOp();

  int getSize() const override;
  int getSizeSimu() const;
  VectorDouble kriging(const VectorDouble& dat) const;
  VectorDouble krigingWithGuess(const VectorDouble& dat, const VectorDouble& guess) const;
  void setMaxIterations(int n) { _solver->setMaxIterations(n); }
  void setTolerance(double tol) { _solver->setTolerance(tol); }
  int getIterations() const { return _solver->getIterations(); }
  double getError() const { return _solver->getError(); }
  VectorDouble computeDriftCoeffs(const VectorDouble& Z,
                                  const MatrixDense& driftMat,
                                  bool verbose = false) const;
  VectorDouble simCond(const VectorDouble& dat) const;
  VectorDouble simNonCond() const;

  const PrecisionOpMulti* getQKriging() const { return _QKriging; }
  const ProjMulti* getProjKriging() const { return _projKriging; }
  const ASimulable* getInvNoise() const { return _invNoise; }
  const PrecisionOpMulti* getQSimu() const { return _QSimu; }
  const ProjMulti* getProjSimu() const { return _projSimu; }

#ifndef SWIG
public:
  int kriging(const constvect inv, vect out) const;
  int krigingWithGuess(const constvect inv,
                       const constvect guess,
                       vect out) const;
  void evalInvCov(const constvect inv, vect result) const;
  void simCond(const constvect data, vect outv) const;
  void simNonCond(vect outv) const;
  virtual double computeLogDetOp(int nbsimu) const;
  double computeQuadratic(const std::vector<double>& x) const;
  double computeTotalLogDet(int nMC) const;

  static int centerDataByDriftMat(VectorDouble& Z,
                                  const MatrixDense& driftMat,
                                  const VectorDouble& driftCoeffs);
  static int centerDataByMeanVec(VectorDouble& Z,
                                 const VectorDouble& meanVec);

protected:
  int _addToDest(const constvect inv, vect outv) const override;

private:
  int _getNDat() const { return _ndat; }
  virtual int _solve(const constvect in, vect out) const;
  int _solveWithGuess(const constvect in,
                      const constvect guess,
                      vect out) const;

  int _buildRhs(const constvect inv) const;
#endif

private:
  void _prepare(bool w1 = true, bool w2 = true) const;

protected:
  const PrecisionOpMulti* const _QKriging;
  const ProjMulti*        const _projKriging;
  const ASimulable*       const _invNoise;
  const PrecisionOpMulti* const _QSimu;
  const ProjMulti*        const _projSimu;
  ALinearOpCGSolver* _solver;

private:
  bool    _noiseToDelete;
  int     _ndat;
  mutable VectorDouble _workdat1; 
  mutable VectorDouble _workdat2;
  mutable VectorDouble _workdat3;
  mutable VectorDouble _workdat4;
  mutable VectorDouble _workNoiseMesh;
  mutable VectorDouble _workNoiseData;
  mutable VectorDouble _rhs;
  mutable VectorDouble _workmesh;
};

/****************************************************************************/

class GSTLEARN_EXPORT SPDEOp : public ASPDEOp,
#ifndef SWIG
  public ALinearOpEigenCG<SPDEOp>
#else
  public virtual ALinearOp
#endif
{
public:
  SPDEOp(const PrecisionOpMulti* const popKriging = nullptr,
         const ProjMulti* const projKriging       = nullptr,
         const ASimulable* const invNoise         = nullptr,
         const PrecisionOpMulti* const popSimu    = nullptr,
         const ProjMulti* const projSimu          = nullptr,
         bool noiseToDelete                       = false)
    : ASPDEOp(popKriging, projKriging, invNoise, popSimu, projSimu, noiseToDelete)
  {
    _solver = new LinearOpCGSolver<SPDEOp>(this);
  }
  virtual ~SPDEOp() = default;


};

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(SPDEOp)
#endif

/****************************************************************************/

#if 0
// To change the algorithm used by SPDEOp, add a new class as below and use
// it instead of SPDEOp:
#ifndef SWIG
#include "LinearOp/ASimulableEigenCG.hpp"
DECLARE_EIGEN_TRAITS(ExampleSPDEOp)
#endif

namespace Eigen {
  namespace internal {
    template<>
    //template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
    /*EIGEN_DONT_INLINE*/ inline void conjugate_gradient(
      const ExampleSPDEOp& /*mat*/,
      const Eigen::Map<const Eigen::VectorXd, 0, Eigen::Stride<0, 0>>& /*rhs*/,
      Eigen::Map<Eigen::VectorXd, 0, Eigen::Stride<0, 0>>& /*x*/,
      const Eigen::IdentityPreconditioner& /*precond*/,
      Index& /*iters*/,
      typename Eigen::Map<Eigen::VectorXd, 0, Eigen::Stride<0, 0>>::RealScalar& /*tol_error*/
    ) {
      messerr("Solver for ExampleSPDEOp");
    }
  }
}

class GSTLEARN_EXPORT ExampleSPDEOp : public ASPDEOp,
#ifndef SWIG
  public ALinearOpEigenCG<ExampleSPDEOp>
#else
  public virtual ALinearOp
#endif
{
public:
  ExampleSPDEOp(const PrecisionOpMulti* const popKriging = nullptr,
                const ProjMulti*        const projKriging = nullptr,
                const ASimulable*       const invNoise = nullptr,
                const PrecisionOpMulti* const popSimu = nullptr,
                const ProjMulti*        const projSimu = nullptr,
                bool  noiseToDelete = false
  ) : ASPDEOp(popKriging, projKriging, invNoise, popSimu, projSimu, noiseToDelete)
  {
    _solver = new LinearOpCGSolver<ExampleSPDEOp>(this);
  }
  virtual ~ExampleSPDEOp() = default;
};

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(ExampleSPDEOp)
#endif

#endif

