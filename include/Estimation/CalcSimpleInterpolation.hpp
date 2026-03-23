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
    CalcSimpleInterpolation(const CalcSimpleInterpolation& r) = delete;
    CalcSimpleInterpolation&
      operator=(const CalcSimpleInterpolation& r) = delete;
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

} // namespace gstlrn
