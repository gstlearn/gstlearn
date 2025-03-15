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

#include "gstlearn_export.hpp"

#include "geoslib_define.h"

#include "Enum/EKrigOpt.hpp"

#include "Calculators/ACalcInterpolator.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

#include "Estimation/CalcKriging.hpp"
class Db;
class DbGrid;
class KrigingSystem;

// TODO : Create KrigingParam which inherits from InterpolatorParam
class GSTLEARN_EXPORT CalcKrigingSimpleCase: public ACalcInterpolator
{
public:
  CalcKrigingSimpleCase(bool flag_est = true, bool flag_std = true, bool flag_varZ = false);
  CalcKrigingSimpleCase(const CalcKrigingSimpleCase &r) = delete;
  CalcKrigingSimpleCase& operator=(const CalcKrigingSimpleCase &r) = delete;
  virtual ~CalcKrigingSimpleCase();

  void setCalcul(const EKrigOpt& calcul);

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  void _storeResultsForExport(const KrigingSystem& ksys);

private:
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;

  EKrigOpt  _calcul;

  VectorString _nameCoord;
  int _iechSingleTarget;
  bool _flagNeighOnly;


  int  _nbNeigh;

  int _iptrEst;
  int _iptrStd;
  int _iptrVarZ;
  int _iptrNeigh;

  Krigtest_Res _ktest;
};

