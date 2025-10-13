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

#include "geoslib_define.h"

#include "Matrix/MatrixDense.hpp"
#include "gstlearn_export.hpp"

#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Model/ModelOptimParam.hpp"

namespace gstlrn
{
class ModelCovList;
class Constraints;
class MatrixDense;
class MatrixSymmetric;

/**
 * \brief
 * Class which, starting from experimental quantities, enables fitting the
 * sills of all Covariance parts of a Model
 */
class GSTLEARN_EXPORT AModelFitSills: public ICloneable
{
public:
  AModelFitSills(ModelCovList* model,
                 const Constraints* constraints = nullptr,
                 const ModelOptimParam& mop     = ModelOptimParam());
  AModelFitSills(const AModelFitSills& m);
  AModelFitSills& operator=(const AModelFitSills& m);
  virtual ~AModelFitSills();

  virtual Id fitSillMatrices() { return 0; }
  Id getNiter() const { return _iterg; }
  void printFitSillSummary(Id niter) const;
  void setTrace(bool trace) { _trace = trace; }
  void setVerbose(bool verbose) { _verbose = verbose; }

protected:
  void _resetInitialSill(std::vector<MatrixSymmetric>& sill) const;
  void _allocateInternalArrays(bool flag_exp = true);
  Id _fitSillMatrices();

private:
  Id _sillFittingIntrinsic();
  Id _goulardWithConstraints();
  Id _goulardWithoutConstraint(Id niter,
                               Id nvar,
                               Id ncova,
                               Id npadir,
                               VectorDouble& wt,
                               VectorDouble& gg,
                               std::vector<MatrixDense>& ge,
                               std::vector<MatrixSymmetric>& sill) const;
  void _storeSillsInModel() const;
  void _optimizeUnderConstraints();
  Id _makeDefinitePositive(Id icov0, double eps = EPSILON12);
  void _initializeGoulard();
  Id _truncateNegativeEigen(Id icov0);
  double _sumSills(Id ivar0, std::vector<MatrixSymmetric>& alpha) const;
  double _calculateScore() const;
  double _calculateScoreMP(Id nvar,
                           Id npadir,
                           VectorDouble& wt,
                           VectorDouble& gg,
                           const MatrixDense& mp) const;
  static Id _combineVariables(Id ivar0, Id jvar0);
  double _minimizeP4(Id icov0,
                     Id ivar0,
                     double xrmax,
                     VectorDouble& xr,
                     std::vector<MatrixSymmetric>& alpha);
  void _updateAlphaDiag(Id icov0,
                        Id ivar0,
                        VectorDouble& xr,
                        std::vector<MatrixSymmetric>& alpha);
  void _updateOtherSills(Id icov0,
                         Id ivar0,
                         std::vector<MatrixSymmetric>& alpha,
                         VectorDouble& xr);
  void _updateCurrentSillGoulard(Id icov0, Id ivar0);
  void _updateCurrentSillDiag(Id icov0,
                              Id ivar0,
                              std::vector<MatrixSymmetric>& alpha,
                              VectorDouble& xr);
  void _updateAlphaNoDiag(Id icov0,
                          Id ivar0,
                          VectorDouble& xr,
                          std::vector<MatrixSymmetric>& alpha);
  bool _convergenceReached(Id iter, double crit, double crit_mem) const;

protected:
  Id _ndim;
  Id _nvar;
  Id _nvs2;
  Id _ncova;
  Id _nbexp;
  Id _npadir;
  VectorDouble _wt;
  VectorDouble _gg;
  VectorDouble _ggc;
  VectorDouble _wtc;
  VectorDouble _wt2;
  VectorDouble _gg2;
  std::vector<VectorDouble> _dd;
  std::vector<MatrixDense> _ge;
  std::vector<MatrixDense> _ge1;
  std::vector<MatrixDense> _ge2;
  std::vector<MatrixSymmetric> _alphau;
  std::vector<MatrixSymmetric> _sill1;
  std::vector<MatrixSymmetric> _sill;

  bool _verbose;
  bool _trace;
  mutable Id _iterg;
  mutable double _score;

  // Storing external pointers or references (not to be deleted)
  ModelCovList* _model;
  const Constraints* _constraints;
  ModelOptimParam _mop;
  CovCalcMode _calcmode;
};
} // namespace gstlrn