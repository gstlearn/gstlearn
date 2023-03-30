/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "geoslib_define.h"

#include "Calculators/ACalcInterpolator.hpp"

class Db;
class DbGrid;
class ANeighParam;

class GSTLEARN_EXPORT CalcSimpleInterpolation: public ACalcInterpolator
{
public:
  CalcSimpleInterpolation();
  CalcSimpleInterpolation(const CalcSimpleInterpolation &r) = delete;
  CalcSimpleInterpolation& operator=(const CalcSimpleInterpolation &r) = delete;
  virtual ~CalcSimpleInterpolation();

  void setFlagMovAve(bool flagMovAve) { _flagMovAve = flagMovAve; }
  void setFlagInvDist(bool flagInvDist) { _flagInvDist = flagInvDist; }
  void setFlagLstSqr(bool flagLstSqr) { _flagLstSqr = flagLstSqr; }

  void setDmax(double dmax) { _dmax = dmax; }
  void setExponent(double exponent) { _exponent = exponent; }
  void setFlagExpand(bool flagExpand) { _flagExpand = flagExpand; }
  void setOrder(int order) { _order = order; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;
  virtual int  _getNVar() const override;

private:
  int    _iattOut;
  bool   _flagMovAve;
  bool   _flagInvDist;
  bool   _flagLstSqr;
  double _exponent;
  bool   _flagExpand;
  double _dmax;
  int    _order;
};

GSTLEARN_EXPORT int inverseDistance(Db *dbin,
                                    Db *dbout,
                                    double exponent = 2.,
                                    bool flag_expand = true,
                                    double dmax = TEST,
                                    const NamingConvention &namconv = NamingConvention(
                                        "InvDist"));
GSTLEARN_EXPORT int movingAverage(Db *dbin,
                                  Db *dbout,
                                  ANeighParam *neighparam,
                                  const NamingConvention &namconv = NamingConvention(
                                      "MovAve"));
GSTLEARN_EXPORT int leastSquares(Db *dbin,
                                 Db *dbout,
                                 ANeighParam *neighparam,
                                 int order = 0,
                                 const NamingConvention &namconv = NamingConvention(
                                     "LstSqr"));

