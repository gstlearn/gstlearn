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

#include "Matrix/MatrixRectangular.hpp"
#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Model/AModelOptim.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"

class Model;
class Constraints;
class MatrixRectangular;
class MatrixSquareSYmmetric;

/**
 * \brief
 * Class which, starting from experimental quantities, enables fitting the
 * sills of all Covariance parts of a Model
 */
class GSTLEARN_EXPORT AModelOptimSills: public AModelOptim
{
public:
  AModelOptimSills(Model* model,
                   Constraints* constraints      = nullptr,
                   const Option_AutoFit& mauto   = Option_AutoFit(),
                   const Option_VarioFit& optvar = Option_VarioFit());
  AModelOptimSills(const AModelOptimSills& m);
  AModelOptimSills& operator=(const AModelOptimSills& m);
  virtual ~AModelOptimSills();

  int fitPerform();

protected:
  void _resetSill(int ncova, std::vector<MatrixSquareSymmetric>& sill) const;
  void _allocateInternalArrays(bool flag_exp = true);

private:
  int _sillFittingIntrinsic(double *crit_arg);
  int _goulardWithConstraints(double *crit_arg);
  int _goulardWithoutConstraint(const Option_AutoFit& mauto,
                                int nvar,
                                int ncova,
                                int npadir,
                                VectorDouble& wt,
                                VectorDouble& gg,
                                std::vector<MatrixRectangular>& ge,
                                std::vector<MatrixSquareSymmetric>& sill,
                                double* crit_arg) const;
  void _storeSillsInModel() const;
  void _optimizeUnderConstraints(double* score);
  int _makeDefinitePositive(int icov0, double eps = EPSILON12);
  void _initializeGoulard();
  int _truncateNegativeEigen(int icov0);
  double _sumSills(int ivar0, std::vector<MatrixSquareSymmetric>& alpha) const;
  double _score();
  static int _combineVariables(int ivar0, int jvar0);
  double _minimizeP4(int icov0,
                     int ivar0,
                     double xrmax,
                     VectorDouble& xr,
                     std::vector<MatrixSquareSymmetric>& alpha);
  void _updateAlphaDiag(int icov0,
                        int ivar0,
                        VectorDouble& xr,
                        std::vector<MatrixSquareSymmetric>& alpha);
  void _updateOtherSills(int icov0,
                         int ivar0,
                         std::vector<MatrixSquareSymmetric>& alpha,
                         VectorDouble& xr);
  void _updateCurrentSillGoulard(int icov0, int ivar0);
  void _updateCurrentSillDiag(int icov0,
                              int ivar0,
                              std::vector<MatrixSquareSymmetric>& alpha,
                              VectorDouble& xr);
  void _updateAlphaNoDiag(int icov0,
                          int ivar0,
                          VectorDouble& xr,
                          std::vector<MatrixSquareSymmetric>& alpha);
  static bool _convergenceReached(const Option_AutoFit& mauto,
                           double crit,
                           double crit_mem);
  void _printResults(double crit) const;

protected:
  int _ndim;
  int _nvar;
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
  std::vector<MatrixRectangular> _ge;
  std::vector<MatrixRectangular> _ge1;
  std::vector<MatrixRectangular> _ge2;
  std::vector<MatrixSquareSymmetric> _alphau;
  std::vector<MatrixSquareSymmetric> _sill1;
  std::vector<MatrixSquareSymmetric> _sill;
};
