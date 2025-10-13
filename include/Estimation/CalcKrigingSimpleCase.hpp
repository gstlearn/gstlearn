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

namespace gstlrn
{
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
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

  void _storeResultsForExport(const KrigingSystemSimpleCase& ksys,
                              KrigingAlgebraSimpleCase& algebra,
                              Id iechout);

private:
  bool _flagEst;
  bool _flagStd;
  bool _flagVarZ;

  EKrigOpt _calcul;

  VectorString _nameCoord;
  Id _iechSingleTarget;

  Id _iptrEst;
  Id _iptrStd;
  Id _iptrVarZ;

  Krigtest_Res _ktest;
};
} // namespace gstlrn