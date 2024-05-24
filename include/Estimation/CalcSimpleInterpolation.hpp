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

#include "Model/Model.hpp"
#include "Calculators/ACalcInterpolator.hpp"

class Db;
class DbGrid;

class GSTLEARN_EXPORT CalcSimpleInterpolation: public ACalcInterpolator
{
public:
  CalcSimpleInterpolation();
  CalcSimpleInterpolation(const CalcSimpleInterpolation &r) = delete;
  CalcSimpleInterpolation& operator=(const CalcSimpleInterpolation &r) = delete;
  virtual ~CalcSimpleInterpolation();

  void setFlagMovAve(bool flagMovAve) { _flagMovAve = flagMovAve; }
  void setFlagMovMed(bool flagMovMed) { _flagMovMed = flagMovMed; }
  void setFlagInvDist(bool flagInvDist) { _flagInvDist = flagInvDist; }
  void setFlagLstSqr(bool flagLstSqr) { _flagLstSqr = flagLstSqr; }
  void setFlagNearest(bool flagNearest) { _flagNearest = flagNearest; }

  void setDmax(double dmax) { _dmax = dmax; }
  void setExponent(double exponent) { _exponent = exponent; }
  void setFlagExpand(bool flagExpand) { _flagExpand = flagExpand; }
  void setOrder(int order) { _order = order; }
  void setFlagEst(bool flagEst) { _flagEst = flagEst; }
  void setFlagStd(bool flagStd) { _flagStd = flagStd; }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;
  virtual bool _postprocess() override;
  virtual void _rollback() override;
  virtual int  _getNVar() const override;

  int _nearest(Db* dbin, Db* dbout, ANeigh* neigh);
  int _movave(Db* dbin, Db* dbout, ANeigh* neigh);
  int _movmed(Db* dbin, Db* dbout, ANeigh* neigh);
  int _lstsqr(Db* dbin, Db* dbout, ANeigh* neigh);
  int _invdist(Db *dbin, Db *dbout);

  void _pointInvdist(Db *dbin, Db *dbout);
  void _gridInvdist(DbGrid *dbin, Db *dbout);

  double _estimCalc(const Db *dbin,
                    const VectorInt &nbgh,
                    const VectorDouble& weights) const;
  double _stdevCalc(Db *dbin,
                    Db *dbout,
                    const VectorInt &nbgh,
                    int iechout,
                    const VectorDouble& weights) const;
  void _saveResults(Db *dbin,
                    Db *dbout,
                    const VectorInt &nbgh,
                    int iech,
                    VectorDouble &weights) const;

private:
  bool   _flagEst;
  bool   _flagStd;
  int    _iattEst;
  int    _iattStd;
  bool   _flagMovAve;
  bool   _flagMovMed;
  bool   _flagInvDist;
  bool   _flagLstSqr;
  bool   _flagNearest;
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
                                    bool flag_est = true,
                                    bool flag_std = false,
                                    Model* model = nullptr,
                                    const NamingConvention &namconv = NamingConvention(
                                        "InvDist"));
GSTLEARN_EXPORT int nearestNeighbor(Db *dbin,
                                    Db *dbout,
                                    bool flag_est = true,
                                    bool flag_std = false,
                                    Model* model = nullptr,
                                    const NamingConvention &namconv = NamingConvention(
                                        "Nearest"));
GSTLEARN_EXPORT int movingAverage(Db *dbin,
                                  Db *dbout,
                                  std::shared_ptr<ANeigh> &neigh,
                                  bool flag_est = true,
                                  bool flag_std = false,
                                  Model *model = nullptr,
                                  const NamingConvention &namconv = NamingConvention(
                                      "MovAve"));
GSTLEARN_EXPORT int movingMedian(Db *dbin,
                                 Db *dbout,
                                 std::shared_ptr<ANeigh>& neigh,
                                 bool flag_est = true,
                                 bool flag_std = false,
                                 Model *model = nullptr,
                                 const NamingConvention &namconv = NamingConvention(
                                     "MovMed"));
GSTLEARN_EXPORT int leastSquares(Db *dbin,
                                 Db *dbout,
                                 std::shared_ptr<ANeigh>& neigh,
                                 int order = 0,
                                 const NamingConvention &namconv = NamingConvention(
                                     "LstSqr"));

