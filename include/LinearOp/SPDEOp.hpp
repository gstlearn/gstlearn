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

#include "LinearOp/ALinearOp.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "LinearOp/ASimulable.hpp"

#ifndef SWIG
#include "LinearOp/ASimulableEigenCG.hpp"
DECLARE_EIGEN_TRAITS(SPDEOp)
#endif

#include "LinearOp/LinearOpCGSolver.hpp"


class PrecisionOpMulti;
class ProjMulti;


class GSTLEARN_EXPORT SPDEOp:
#ifndef SWIG
  public ASimulableEigenCG<SPDEOp>
#else
  public ASimulable
#endif
{

public:
  SPDEOp(const PrecisionOpMulti* const pop      = nullptr, 
         const ProjMulti*        const proj     = nullptr,
         const ASimulable*       const invNoise = nullptr,
         bool  noiseToDelete                    = false);
  virtual ~SPDEOp();

  int getSize() const override;
  VectorDouble kriging(const VectorDouble& dat) const;
  VectorDouble krigingWithGuess(const VectorDouble& dat, const VectorDouble& guess) const;
  void setMaxIterations(int n) {_solver.setMaxIterations(n);}
  void setTolerance(double tol) {_solver.setTolerance(tol);}
  int  getIterations() const { return _solver.getIterations();}
  double getError() const { return  _solver.getError();}
  VectorDouble computeDriftCoeffs(const VectorDouble& Z, 
                                  const MatrixRectangular& drifts) const;
#ifndef SWIG
public:
  int kriging(const constvect& inv,
                    vect& out) const;
  int krigingWithGuess(const constvect& inv,
                       const constvect& guess,
                             vect& out) const;
  void evalInvCov(const constvect &inv, vect& result) const;
 
protected:
  int _addToDest(const constvect& inv, vect& outv) const override;
  int _addSimulateToDest(const constvect& whitenoise, vect& outv) const override;

private: 
  int _getNDat() const {return _ndat;}
  virtual int _solve(const constvect& in,vect& out) const;
  int _solveWithGuess(const constvect& in,const constvect &guess,vect& out) const;

  int _buildRhs(const constvect& inv) const;
#endif

private:
  void _prepare(bool w1 = true, bool w2 = true) const;
#ifndef SWIG
private:
  virtual int _addToDestImpl(const constvect& inv, vect& outv) const;
#endif

protected:
  const PrecisionOpMulti* const _Q;
  const ProjMulti*        const _Proj;
  const ASimulable*       const _invNoise;

private:
  bool    _noiseToDelete;
  int     _ndat;
  mutable VectorDouble _workdat1; 
  mutable VectorDouble _workdat2;
  mutable VectorDouble _rhs;
  mutable VectorDouble _workmesh;
  mutable LinearOpCGSolver<SPDEOp> _solver;

};

#ifndef SWIG
DECLARE_EIGEN_PRODUCT(SPDEOp)
#endif