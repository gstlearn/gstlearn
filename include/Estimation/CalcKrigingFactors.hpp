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

#include "Calculators/ACalcInterpolator.hpp"

class Db;
class DbGrid;
class KrigingSystem;

class GSTLEARN_EXPORT CalcKrigingFactors: public ACalcInterpolator
{
public:
  CalcKrigingFactors(bool flag_est = true, bool flag_std = true);
  CalcKrigingFactors(const CalcKrigingFactors &r) = delete;
  CalcKrigingFactors& operator=(const CalcKrigingFactors &r) = delete;
  virtual ~CalcKrigingFactors();

  void setIuidFactors(const VectorInt& iuidFactors) { _iuidFactors = iuidFactors; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  int  _getNFactors() const;
  void _storeResultsForExport(const KrigingSystem& ksys);
  bool _hasChangeSupport() const;

private:
  bool _flagEst;
  bool _flagStd;

  VectorString _nameCoord;

  int _iptrEst;
  int _iptrStd;

  VectorInt _iuidFactors;

  Model* _modelLocal;
};

GSTLEARN_EXPORT int krigingFactors(Db* dbin,
                                   Db* dbout,
                                   Model* model,
                                   ANeigh* neigh,
                                   bool flag_est                   = true,
                                   bool flag_std                   = true,
                                   const KrigOpt& krigopt          = KrigOpt(),
                                   const NamingConvention& namconv = NamingConvention("KD"));
