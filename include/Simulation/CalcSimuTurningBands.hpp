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

#include "Covariances/CovAniso.hpp"
#include "Enum/ECov.hpp"
#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Model/Model.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/TurningBandDirection.hpp"

#include "geoslib_define.h"

namespace gstlrn
{
  class Model;
  class ANeigh;
  class TurningBandOperate;

  class GSTLEARN_EXPORT CalcSimuTurningBands: public ACalcSimulation
  {
  public:
    CalcSimuTurningBands(
      Id nbsimu = 0,
      Id nbtuba = 0,
      bool flag_check = false,
      Id seed = 4324324);
    CalcSimuTurningBands(const CalcSimuTurningBands& r) = delete;
    CalcSimuTurningBands& operator=(const CalcSimuTurningBands& r) = delete;
    virtual ~CalcSimuTurningBands();

    Id getNBtuba() const { return _nbtuba; }

    void setNBtuba(Id nbtuba) { _nbtuba = nbtuba; }

    Id getNDirs() const { return static_cast<Id>(_codirs.size()); }

    Id simulate(
      Db* dbin,
      Db* dbout,
      Model* model,
      ANeigh* neigh,
      Id icase,
      Id flag_bayes = false,
      bool flag_pgs = false,
      bool flag_gibbs = false,
      bool flag_dgm = false);
    Id simulatePotential(
      Db* dbiso,
      Db* dbgrd,
      Db* dbtgt,
      Db* dbout,
      ModelGeneric* model,
      double delta);

    bool isFlagCheck() const { return _flagCheck; }

    void setFlagCheck(bool flag_check) { _flagCheck = flag_check; }

    void setIcase(Id icase) { _icase = icase; }

    Id getNbtuba() const { return _nbtuba; }

    void setNbtuba(Id nbtuba) { _nbtuba = nbtuba; }

    void setFlagAllocationAlreadyDone(Id flag)
    {
      _flagAllocationAlreadyDone = flag;
    }

  private:
    bool _check() override;
    bool _preprocess() override;
    bool _run() override;
    bool _postprocess() override;
    void _rollback() override;
    bool _simulate();
    Id _computeTB(
      Db* db,
      Id isimu,
      const VectorBool& activeArray,
      VectorVectorDouble& tab) override;
    bool _resizeTB();
    void _computePoint(
      Db* db,
      const CovAniso* cova,
      const ECov& type,
      Id isimu,
      Id is,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);
    void _computeGrid(
      DbGrid* db,
      const CovAniso* cova,
      const ECov& type,
      Id isimu,
      Id is,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);
    void _computeNugget(
      Db* db,
      Id isimu,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);

    // Turning bands specific methods
    Id _getIBS(Id isimu, Id is, Id ib) const;

    Id _getIcase() const override { return _icase; }

    void _setCodirAng(Id ibs, Id idir, double value)
    {
      _codirs[ibs].setAng(idir, value);
    }

    void _setCodirTmin(Id ibs, double value) { _codirs[ibs].setTmin(value); }

    void _setCodirTmax(Id ibs, double value) { _codirs[ibs].setTmax(value); }

    void _setCodirScale(Id ibs, double value) { _codirs[ibs].setScale(value); }

    void _setCodirT00(Id ibs, double value) { _codirs[ibs].setT00(value); }

    void _setCodirDXP(Id ibs, double value) { _codirs[ibs].setDXP(value); }

    void _setCodirDYP(Id ibs, double value) { _codirs[ibs].setDYP(value); }

    void _setCodirDZP(Id ibs, double value) { _codirs[ibs].setDZP(value); }

    VectorDouble _getCodirAng(Id ibs) const { return _codirs[ibs].getAng(); }

    double _getCodirAng(Id ibs, Id idir) const
    {
      return _codirs[ibs].getAng(idir);
    }

    double _getCodirScale(Id ibs) { return _codirs[ibs].getScale(); }

    double _getCodirT00(Id ibs) const { return _codirs[ibs].getT00(); }

    double _getCodirDXP(Id ibs) const { return _codirs[ibs].getDXP(); }

    double _getCodirDYP(Id ibs) const { return _codirs[ibs].getDYP(); }

    double _getCodirDZP(Id ibs) const { return _codirs[ibs].getDZP(); }

    double _getCodirTmin(Id ibs) const { return _codirs[ibs].getTmin(); }

    double _getCodirTmax(Id ibs) const { return _codirs[ibs].getTmax(); }

    Id _getAddressBand(Id ivar, Id is, Id ib, Id isimu) const;
    void _setSeedBand(Id ivar, Id is, Id ib, Id isimu, Id seed);
    Id _getSeedBand(Id ivar, Id is, Id ib, Id isimu) const;

    void _rotateDirections(double a[3], double theta);
    Id _generateDirections(const Db* dbout);
    void _minmax(const Db* db);
    void _setDensity();
    static ECov _particularCase(const CovAniso* cova, double eps = EPSILON7);
    Id _initializeSeedBands();
    void _normalizeForBands(
      const Db* db,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);
    Id _getCorrec(
      const ECov& type,
      Id is,
      Id ibs,
      TurningBandOperate& operTB,
      double& correc);
    static double _getScale(double alpha, double scale);
    static double _getScaleKB(double param, double scale);
    void _migrationInit(
      Id ibs,
      Id is,
      double scale,
      TurningBandOperate& operTB,
      double eps = EPSILON5);
    double _dilutionInit(Id ibs, Id is, TurningBandOperate& operTB);
    double _spectralInit(Id ibs, Id is, TurningBandOperate& operTB);
    double _power1DInit(Id ibs, Id is, TurningBandOperate& operTB);
    double _spline1DInit(Id ibs, Id k, TurningBandOperate& operTB);
    double _irfProcessInit(Id ibs, Id is, TurningBandOperate& operTB);

    static double _irfCorrec(const ECov& type, double theta1, double scale);
    void _getOmegaPhi(
      Id ibs,
      TurningBandOperate& operTB,
      double* cxp,
      double* sxp,
      double* cyp,
      double* syp,
      double* czp,
      double* szp,
      double* c0z,
      double* s0z);
    void _spreadRegularOnGrid(
      const DbGrid* dbgrid,
      const CovAniso* cova,
      Id ibs,
      double correc,
      TurningBandOperate& operTB,
      const VectorBool& activeArray,
      VectorDouble& tab);
    void _spreadRegularOnPoint(
      const Db* db,
      const CovAniso* cova,
      Id ibs,
      double correc,
      TurningBandOperate& operTB,
      const VectorBool& activeArray,
      VectorDouble& tab);
    void _spreadSpectralOnGrid(
      const DbGrid* dbgrid,
      const CovAniso* cova,
      Id ibs,
      double correc,
      TurningBandOperate& operTB,
      const VectorBool& activeArray,
      VectorDouble& tab);
    void _spreadSpectralOnPoint(
      const Db* db,
      const CovAniso* cova,
      Id ibs,
      double correc,
      TurningBandOperate& operTB,
      const VectorBool& activeArray,
      VectorDouble& tab);

    // Debugging methods
    void _dumpBands() const;
    void _dumpSeeds() const;

  private:
    Id _nbtuba;
    Id _iattOut;
    Id _icase;
    bool _flagCheck;
    bool _flagAllocationAlreadyDone;
    VectorString _nameCoord;
    Id _npointSimulated;
    double _field;
    double _theta;
    VectorInt _seedBands;
    std::vector<TurningBandDirection> _codirs;
    Model*
      _modelLocal; // Conversion of getModel() into a Model (more than ModelGeneric)
  };

} // namespace gstlrn
