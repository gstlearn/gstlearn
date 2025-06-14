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

#include "Basic/VectorNumT.hpp"
#include "Model/ModelOptimParam.hpp"
#include "Basic/ICloneable.hpp"
#include "Covariances/CovCalcMode.hpp"

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

  virtual int fitSills(bool verbose = false, bool trace = false) { DECLARE_UNUSED(verbose, trace); return 0; }
  int getNiter() const { return _iterg; }
  void printFitSillSummary(int niter) const;

protected:
  void _resetInitialSill(std::vector<MatrixSymmetric>& sill) const;
  void _allocateInternalArrays(bool flag_exp = true);
  int  _fitSills(bool verbose = false, bool trace = false);

private:
  int _sillFittingIntrinsic();
  int _goulardWithConstraints();
  int _goulardWithoutConstraint(int niter,
                                int nvar,
                                int ncova,
                                int npadir,
                                VectorDouble& wt,
                                VectorDouble& gg,
                                std::vector<MatrixDense>& ge,
                                std::vector<MatrixSymmetric>& sill) const;
  void _storeSillsInModel() const;
  void _optimizeUnderConstraints();
  int _makeDefinitePositive(int icov0, double eps = EPSILON12);
  void _initializeGoulard();
  int _truncateNegativeEigen(int icov0);
  double _sumSills(int ivar0, std::vector<MatrixSymmetric>& alpha) const;
  double _score();
  static int _combineVariables(int ivar0, int jvar0);
  double _minimizeP4(int icov0,
                     int ivar0,
                     double xrmax,
                     VectorDouble& xr,
                     std::vector<MatrixSymmetric>& alpha);
  void _updateAlphaDiag(int icov0,
                        int ivar0,
                        VectorDouble& xr,
                        std::vector<MatrixSymmetric>& alpha);
  void _updateOtherSills(int icov0,
                         int ivar0,
                         std::vector<MatrixSymmetric>& alpha,
                         VectorDouble& xr);
  void _updateCurrentSillGoulard(int icov0, int ivar0);
  void _updateCurrentSillDiag(int icov0,
                              int ivar0,
                              std::vector<MatrixSymmetric>& alpha,
                              VectorDouble& xr);
  void _updateAlphaNoDiag(int icov0,
                          int ivar0,
                          VectorDouble& xr,
                          std::vector<MatrixSymmetric>& alpha);
  bool _convergenceReached(int iter, double crit, double crit_mem) const;

protected:
  int _ndim;
  int _nvar;
  int _nvs2;
  int _ncova;
  int _nbexp;
  int _npadir;
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
  mutable int _iterg;
  mutable double _crit;

  // Storing external pointers or references (not to be deleted)
  ModelCovList*      _model;
  const Constraints* _constraints;
  ModelOptimParam    _mop;
  CovCalcMode        _calcmode;
};
