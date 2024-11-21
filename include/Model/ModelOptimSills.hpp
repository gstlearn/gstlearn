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
#include "Model/ModelOptimVario.hpp"

class Model;
class Vario;
class Constraints;
class Option_AutoFit;
class Option_VarioFit;
class MatrixRectangular;
class MatrixSquareSYmmetric;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * sills of all Covariance parts of a Model
 */
class GSTLEARN_EXPORT ModelOptimSills: public ModelOptimVario
{
public:
  ModelOptimSills();
  ModelOptimSills(const ModelOptimSills& m);
  ModelOptimSills& operator=(const ModelOptimSills& m);
  virtual ~ModelOptimSills();

  int fit(Vario* vario,
          Model* model,
          Constraints* constraints,
          const Option_AutoFit* mauto,
          const Option_VarioFit* optvar);

private:
  int  _getDimensions();
  void _allocateInternalArrays(bool flag_exp = true);
  int  _constraintsCheck();
  void _computeGg();
  void _computeGe();
  void _compressArray(const VectorDouble& tabin, VectorDouble& tabout);
  void _preparGoulardVario();
  int _goulardWithoutConstraint(const Option_AutoFit* mauto,
                                int nvar,
                                int ncova,
                                int npadir,
                                VectorDouble& wt,
                                VectorDouble& gg,
                                std::vector<MatrixRectangular>& ge,
                                std::vector<MatrixSquareSymmetric>& sill,
                                double* crit_arg) const;
  int _goulardWithConstraints();
  void _initializeGoulard();
  int _makeDefinitePositive(int icov0, double eps = EPSILON12);
  int _optimizeUnderConstraints(double* score);
  int _truncateNegativeEigen(int icov0);
  double _score();
  double _sumSills(int ivar0, std::vector<MatrixSquareSymmetric>& alpha) const;
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
  static int _combineVariables(int ivar0, int jvar0);
  void _resetSill(int ncova, std::vector<MatrixSquareSymmetric>& sill) const;
  void _storeSillsInModel() const;
  int _sillFittingIntrinsic();

private:
  Constraints* _constraints;
  const Option_AutoFit* _mauto;
  const Option_VarioFit* _optvar;

  int _ndim;
  int _nvar;
  int _ncova;
  int _nbexp;
  int _npadir;
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
