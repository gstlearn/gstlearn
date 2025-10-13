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
#include "LinearOp/ALinearOpEigenCG.hpp"
#include "LinearOp/ASimulable.hpp"
#include "LinearOp/LinearOpCGSolver.hpp"
#include "Matrix/MatrixDense.hpp"

#ifndef SWIG
DECLARE_EIGEN_TRAITS(SPDEOp)
#endif

namespace gstlrn
{

class ProjMulti;
class ALinearOp;
class PrecisionOpMulti;
class Chebychev;

class GSTLEARN_EXPORT ASPDEOp: public virtual ALinearOp
{
public:
  ASPDEOp(const PrecisionOpMulti* const popKriging = nullptr,
          const ProjMulti* const projInKriging     = nullptr,
          const ASimulable* invNoise               = nullptr,
          const PrecisionOpMulti* const popSimu    = nullptr,
          const ProjMulti* const projInSimu        = nullptr);
  virtual ~ASPDEOp();

  virtual VectorDouble stdev(const VectorDouble& dat,
                             Id nMC                 = 1,
                             Id seed                = 134343,
                             const ProjMulti* projK = nullptr,
                             const ProjMulti* projS = nullptr) const;

  Id getSize() const override;
  Id getSizeSimu() const;
  Id getIterations() const { return _solver->getIterations(); }
  double getError() const { return _solver->getError(); }

  void setMaxIterations(Id n) { _solver->setMaxIterations(n); }
  void setTolerance(double tol) { _solver->setTolerance(tol); }

  VectorDouble kriging(const VectorDouble& dat, const ProjMulti* proj = nullptr) const;
  VectorDouble krigingWithGuess(const VectorDouble& dat, const VectorDouble& guess) const;

  VectorDouble computeDriftCoeffs(const VectorDouble& Z,
                                  const MatrixDense& driftMat,
                                  bool verbose = false) const;
  VectorDouble simCond(const VectorDouble& dat,
                       const ProjMulti* projK = nullptr,
                       const ProjMulti* projS = nullptr) const;
  VectorDouble simNonCond(const ProjMulti* proj = nullptr) const;

  const PrecisionOpMulti* getQKriging() const { return _QKriging; }
  const ProjMulti* getProjKriging() const { return _projInKriging; }
  const ASimulable* getInvNoise() const { return _invNoise; }
  const PrecisionOpMulti* getQSimu() const { return _QSimu; }
  const ProjMulti* getProjInSimu() const { return _projInSimu; }

#ifndef SWIG

public:
  Id krigingWithGuess(const constvect inv,
                      const constvect guess,
                      vect out) const;
  void evalInvCov(const constvect inv, vect result) const;
  void simCond(const constvect data, vect outv) const;
  void simNonCond(vect outv) const;
#endif
  virtual double computeLogDetOp(Id nbsimu = 1) const;
  double computeQuadratic(const std::vector<double>& x) const;
  double computeTotalLogDet(Id nMC = 5, Id seed = 13132) const;
  double computeLogDetQ(Id nMC = 5) const;
  double computeLogDetInvNoise() const;
  static Id centerDataByDriftMat(VectorDouble& Z,
                                 const MatrixDense& driftMat,
                                 const VectorDouble& driftCoeffs);
  static Id centerDataByMeanVec(VectorDouble& Z,
                                const VectorDouble& meanVec);
  void setVerbose(bool v) { _verbose = v; }

#ifndef SWIG

protected:
  Id _addToDest(const constvect inv, vect outv) const override;

private:
  std::pair<double, double> _computeRangeEigenVal() const;
  void _preparePoly(Chebychev& logPoly) const;
  Id _kriging(const constvect inv, vect out) const;
  void _simNonCond(vect outv) const;
  void _simCond(const constvect data, vect outvK, vect outvS) const;
  Id _getNDat() const { return _ndat; }
  virtual Id _solve(const constvect in, vect out) const;
  Id _solveWithGuess(const constvect in,
                     const constvect guess,
                     vect out) const;
  Id _buildRhs(const constvect inv) const;
#endif

private:
  void _prepare(bool w1 = true, bool w2 = true) const;

protected:
  const PrecisionOpMulti* const _QKriging;
  const ProjMulti* const _projInKriging;
  const ASimulable* const _invNoise;
  const PrecisionOpMulti* const _QSimu;
  const ProjMulti* const _projInSimu;
  ALinearOpCGSolver* _solver;
  bool _verbose;

private:
  Id _ndat;
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

class GSTLEARN_EXPORT SPDEOp: public ASPDEOp,
#ifndef SWIG
                              public ALinearOpEigenCG<SPDEOp>
#else
                              public virtual ALinearOp
#endif
{
public:
  SPDEOp(const PrecisionOpMulti* const popKriging = nullptr,
         const ProjMulti* const projInKriging     = nullptr,
         const ASimulable* invNoise               = nullptr,
         const PrecisionOpMulti* const popSimu    = nullptr,
         const ProjMulti* const projInSimu        = nullptr)
    : ASPDEOp(popKriging, projInKriging, invNoise, popSimu, projInSimu)
  {
    _solver = new LinearOpCGSolver<SPDEOp>(this);
  }
  virtual ~SPDEOp() = default;
};

} // namespace gstlrn

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(SPDEOp)
#endif

/****************************************************************************/

#if 0
// To change the algorithm used by SPDEOp, add a new class as below and use
// it instead of SPDEOp:
#  ifndef SWIG
#    include "LinearOp/ASimulableEigenCG.hpp"
DECLARE_EIGEN_TRAITS(ExampleSPDEOp)
#  endif

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
#  ifndef SWIG
  public ALinearOpEigenCG<ExampleSPDEOp>
#  else
  public virtual ALinearOp
#  endif
{
public:
  ExampleSPDEOp(const PrecisionOpMulti* const popKriging = nullptr,
                const ProjMulti*        const projInKriging = nullptr,
                const ASimulable*       const invNoise = nullptr,
                const PrecisionOpMulti* const popSimu = nullptr,
                const ProjMulti*        const projInSimu = nullptr)
  ) : ASPDEOp(popKriging, projInKriging, invNoise, popSimu, projInSimu)
  {
    _solver = new LinearOpCGSolver<ExampleSPDEOp>(this);
  }
  virtual ~ExampleSPDEOp() = default;
};

#  ifndef SWIG
DECLARE_EIGEN_PRODUCT(ExampleSPDEOp)
#  endif

#endif
