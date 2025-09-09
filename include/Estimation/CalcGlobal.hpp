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

#include "Model/ModelGeneric.hpp"
#include "gstlearn_export.hpp"

#include "geoslib_define.h"

#include "Calculators/ACalcInterpolator.hpp"

namespace gstlrn
{

class Db;
class DbGrid;
class KrigingSystem;
class GSTLEARN_EXPORT Global_Result
{
public:
  Id ntot;              // Total Number of Data
  Id np;                // Number of active Data
  Id ng;                // Number of grid nodes for Domain discretization
  double surface;       // Surface of Domain
  double zest;          // Estimate
  double sse;           // Standard deviation of estimation
  double cvgeo;         // Coefficient of Variation
  double cvv;           // Variance of Domain
  VectorDouble weights; // Weights attached to data

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;
};

class GSTLEARN_EXPORT CalcGlobal: public ACalcInterpolator
{
public:
  CalcGlobal(Id ivar0     = 0,
             bool verbose = false);
  CalcGlobal(const CalcGlobal& r)            = delete;
  CalcGlobal& operator=(const CalcGlobal& r) = delete;
  virtual ~CalcGlobal();

  void setFlagArithmetic(bool flagArithmetic) { _flagArithmetic = flagArithmetic; }
  void setFlagKriging(bool flagKriging) { _flagKriging = flagKriging; }

  Global_Result getGRes() const { return _gRes; }

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

  Id _globalKriging();
  Id _globalArithmetic();

private:
  bool _flagArithmetic;
  bool _flagKriging;
  Id _ivar0;
  bool _verbose;
  Model* _modelLocal;

  Global_Result _gRes;
};

GSTLEARN_EXPORT Global_Result global_arithmetic(Db* dbin,
                                                DbGrid* dbgrid,
                                                ModelGeneric* model,
                                                Id ivar0     = 0,
                                                bool verbose = false);
GSTLEARN_EXPORT Global_Result global_kriging(Db* dbin,
                                             Db* dbout,
                                             ModelGeneric* model,
                                             Id ivar0     = 0,
                                             bool verbose = false);
} // namespace gstlrn