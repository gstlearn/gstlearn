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

#include "Calculators/ACalcInterpolator.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{

class Db;
class DbGrid;
class KrigingSystem;

class GSTLEARN_EXPORT CalcKrigingGradient: public ACalcInterpolator
{
public:
  CalcKrigingGradient(bool flag_est      = true,
                      bool flag_std      = true,
                      double ball_radius = 0.);
  CalcKrigingGradient(const CalcKrigingGradient& r)            = delete;
  CalcKrigingGradient& operator=(const CalcKrigingGradient& r) = delete;
  virtual ~CalcKrigingGradient();

  void setBallRadius(double ball_radius) { _ballRadius = ball_radius; }
  void setFlagForceNumeric(bool status) { _flagForceNumeric = status; }

private:
  void _updateDbin();
  Id _updateModel();

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

private:
  Db* _dbGradient;
  ModelGeneric* _modelGradient;
  double _ballRadius;
  bool _flagEst;
  bool _flagStd;
  bool _flagForceNumeric;
};

GSTLEARN_EXPORT Id krigingGradient(Db* dbin,
                                   Db* dbout,
                                   ModelGeneric* model,
                                   ANeigh* neigh,
                                   bool flag_est                   = true,
                                   bool flag_std                   = true,
                                   double ball_radius              = 0.01,
                                   bool flagForceNumeric           = false,
                                   const NamingConvention& namconv = NamingConvention("KrigGradient"));

} // namespace gstlrn
