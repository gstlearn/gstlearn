/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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

  void setCalcul(const EKrigOpt &calcul) { _calcul = calcul; }
  void setNdisc(const VectorInt &ndisc) { _ndisc = ndisc; }
  void setIuidFactors(const VectorInt &iuidFactors) { _iuidFactors = iuidFactors; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;

  int _getNFactors() const;
  void _storeResultsForExport(const KrigingSystem& ksys);

private:
  bool _flagEst;
  bool _flagStd;

  EKrigOpt  _calcul;
  VectorInt _ndisc;

  int _iptrEst;
  int _iptrStd;

  VectorInt _iuidFactors;
};

GSTLEARN_EXPORT int KrigingFactors(Db *dbin,
                                   Db *dbout,
                                   Model *model,
                                   ANeighParam *neighparam,
                                   const EKrigOpt &calcul = EKrigOpt::fromKey(
                                       "PONCTUAL"),
                                   const VectorInt &ndisc = VectorInt(),
                                   bool flag_est = true,
                                   bool flag_std = true,
                                   const NamingConvention &namconv = NamingConvention(
                                       "KD"));
