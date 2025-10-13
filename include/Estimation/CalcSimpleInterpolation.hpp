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

namespace gstlrn
{

class Db;
class DbGrid;
class Model;

class GSTLEARN_EXPORT CalcSimpleInterpolation: public ACalcInterpolator
{
public:
  CalcSimpleInterpolation();
  CalcSimpleInterpolation(const CalcSimpleInterpolation& r)            = delete;
  CalcSimpleInterpolation& operator=(const CalcSimpleInterpolation& r) = delete;
  virtual ~CalcSimpleInterpolation();

  void setFlagMovAve(bool flagMovAve) { _flagMovAve = flagMovAve; }
  void setFlagMovMed(bool flagMovMed) { _flagMovMed = flagMovMed; }
  void setFlagInvDist(bool flagInvDist) { _flagInvDist = flagInvDist; }
  void setFlagLstSqr(bool flagLstSqr) { _flagLstSqr = flagLstSqr; }
  void setFlagNearest(bool flagNearest) { _flagNearest = flagNearest; }

  void setDmax(double dmax) { _dmax = dmax; }
  void setExponent(double exponent) { _exponent = exponent; }
  void setFlagExpand(bool flagExpand) { _flagExpand = flagExpand; }
  void setOrder(Id order) { _order = order; }
  void setFlagEst(bool flagEst) { _flagEst = flagEst; }
  void setFlagStd(bool flagStd) { _flagStd = flagStd; }

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;
  bool _postprocess() override;
  void _rollback() override;

  Id _nearest(Db* dbin, Db* dbout, ANeigh* neigh);
  Id _movave(Db* dbin, Db* dbout, ANeigh* neigh);
  Id _movmed(Db* dbin, Db* dbout, ANeigh* neigh);
  Id _lstsqr(Db* dbin, Db* dbout, ANeigh* neigh) const;
  Id _invdist(Db* dbin, Db* dbout);

  void _pointInvdist(Db* dbin, Db* dbout);
  void _gridInvdist(DbGrid* dbin, Db* dbout);

  static double _estimCalc(const Db* dbin,
                           const VectorInt& nbgh,
                           const VectorDouble& weights);
  double _stdevCalc(Db* dbin,
                    Db* dbout,
                    const VectorInt& nbgh,
                    Id iechout,
                    const VectorDouble& weights) const;
  void _saveResults(Db* dbin,
                    Db* dbout,
                    const VectorInt& nbgh,
                    Id iech,
                    VectorDouble& weights) const;

private:
  bool _flagEst;
  bool _flagStd;
  Id _iattEst;
  Id _iattStd;
  bool _flagMovAve;
  bool _flagMovMed;
  bool _flagInvDist;
  bool _flagLstSqr;
  bool _flagNearest;
  double _exponent;
  bool _flagExpand;
  double _dmax;
  Id _order;
};

GSTLEARN_EXPORT Id inverseDistance(Db* dbin,
                                    Db* dbout,
                                    double exponent                 = 2.,
                                    bool flag_expand                = true,
                                    double dmax                     = TEST,
                                    bool flag_est                   = true,
                                    bool flag_std                   = false,
                                    Model* model                    = nullptr,
                                    const NamingConvention& namconv = NamingConvention(
                                      "InvDist"));
GSTLEARN_EXPORT Id nearestNeighbor(Db* dbin,
                                    Db* dbout,
                                    bool flag_est                   = true,
                                    bool flag_std                   = false,
                                    Model* model                    = nullptr,
                                    const NamingConvention& namconv = NamingConvention(
                                      "Nearest"));
GSTLEARN_EXPORT Id movingAverage(Db* dbin,
                                  Db* dbout,
                                  ANeigh* neigh,
                                  bool flag_est                   = true,
                                  bool flag_std                   = false,
                                  Model* model                    = nullptr,
                                  const NamingConvention& namconv = NamingConvention(
                                    "MovAve"));
GSTLEARN_EXPORT Id movingMedian(Db* dbin,
                                 Db* dbout,
                                 ANeigh* neigh,
                                 bool flag_est                   = true,
                                 bool flag_std                   = false,
                                 Model* model                    = nullptr,
                                 const NamingConvention& namconv = NamingConvention(
                                   "MovMed"));
GSTLEARN_EXPORT Id leastSquares(Db* dbin,
                                 Db* dbout,
                                 ANeigh* neigh,
                                 Id order                       = 0,
                                 const NamingConvention& namconv = NamingConvention(
                                   "LstSqr"));
} // namespace gstlrn