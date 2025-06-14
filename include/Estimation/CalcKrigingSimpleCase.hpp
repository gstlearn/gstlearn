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

#include "Estimation/KrigingAlgebraSimpleCase.hpp"
#include "gstlearn_export.hpp"

#include "Enum/EKrigOpt.hpp"

#include "Calculators/ACalcInterpolator.hpp"

#include "Estimation/CalcKriging.hpp"
class Db;
class DbGrid;
class KrigingSystemSimpleCase;

// TODO : Create KrigingParam which inherits from InterpolatorParam
class GSTLEARN_EXPORT CalcKrigingSimpleCase: public ACalcInterpolator
{
public:
  CalcKrigingSimpleCase(bool flag_est = true, bool flag_std = true, bool flag_varZ = false);
  CalcKrigingSimpleCase(const CalcKrigingSimpleCase& r)            = delete;
  CalcKrigingSimpleCase& operator=(const CalcKrigingSimpleCase& r) = delete;
  virtual ~CalcKrigingSimpleCase();

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  void _storeResultsForExport(const KrigingSystemSimpleCase& ksys,
                              KrigingAlgebraSimpleCase& algebra,
                              int iechout);

private:
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;

  EKrigOpt _calcul;

  VectorString _nameCoord;
  int _iechSingleTarget;

  int _iptrEst;
  int _iptrStd;
  int _iptrVarZ;

  Krigtest_Res _ktest;
};
