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

class Db;
class DbGrid;
class KrigingSystem;

class GSTLEARN_EXPORT Global_Result
{
public:
  int ntot; // Total Number of Data
  int np;   // Number of active Data
  int ng;   // Number of grid nodes for Domain discretization
  double surface; // Surface of Domain
  double zest;    // Estimate
  double sse;     // Standard deviation of estimation
  double cvgeo;   // Coefficient of Variation
  double cvv;     // Variance of Domain
  VectorDouble weights; // Weights attached to data

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;
};

class GSTLEARN_EXPORT CalcGlobal: public ACalcInterpolator
{
public:
  CalcGlobal(int ivar0 = 0,
             bool verbose = false);
  CalcGlobal(const CalcGlobal &r) = delete;
  CalcGlobal& operator=(const CalcGlobal &r) = delete;
  virtual ~CalcGlobal();

  void setFlagArithmetic(bool flagArithmetic) { _flagArithmetic = flagArithmetic; }
  void setFlagKriging(bool flagKriging) { _flagKriging = flagKriging; }

  Global_Result getGRes() const { return _gRes; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  int _globalKriging();
  int _globalArithmetic();

private:
  bool _flagArithmetic;
  bool _flagKriging;
  int    _ivar0;
  bool   _verbose;
  Model* _modelLocal;

  Global_Result _gRes;
};

GSTLEARN_EXPORT Global_Result global_arithmetic(Db *dbin,
                                                DbGrid *dbgrid,
                                                ModelGeneric *model,
                                                int ivar0,
                                                bool verbose);
GSTLEARN_EXPORT Global_Result global_kriging(Db *dbin,
                                             Db *dbout,
                                             ModelGeneric *model,
                                             int ivar0,
                                             bool verbose);
