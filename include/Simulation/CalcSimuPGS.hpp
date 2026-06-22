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

#include "Simulation/ACalcSimulation.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include <array>

namespace gstlrn
{
  class RuleProp;
  class Model;
  class ANeigh;
  class Rule;

  class GSTLEARN_EXPORT CalcSimuPGS: public ACalcSimulation
  {
  public:
    CalcSimuPGS(
      Id nbsimu = 0,
      Id nbtuba = 100,
      Id gibbsNBurn = 10,
      Id gibbsNIter = 100,
      double percent = 5.0,
      Id seed = 43234,
      bool verbose = false);
    CalcSimuPGS(const CalcSimuPGS& r) = delete;
    CalcSimuPGS& operator=(const CalcSimuPGS& r) = delete;
    virtual ~CalcSimuPGS() = default;

    void setFlagGaus(bool flag) { _flagGaus = flag; }

    void setFlagProp(bool flag) { _flagProp = flag; }

    void setFlagCheck(bool flag) { _flagCheck = flag; }

    void setFlagShow(bool flag) { _flagShow = flag; }

    void setNbTuba(Id nbtuba) { _nbtuba = nbtuba; }

    void setGibbsNBurn(Id gibbsNBurn) { _gibbsNBurn = gibbsNBurn; }

    void setGibbsNIter(Id gibbsNIter) { _gibbsNIter = gibbsNIter; }

    void setPercent(double percent) { _percent = percent; }

    void setRuleProp(RuleProp* ruleprop) { _ruleprop = ruleprop; }

    void setModel(Id ipgs, Id igrf, Model* model)
    {
      if (ipgs < 2 && igrf < 2) _models[ipgs][igrf] = model;
    }

  private:
    bool _check() override;
    bool _preprocess() override;
    bool _run() override;
    bool _postprocess() override;
    void _rollback() override;

    bool _isEnvironmentValid(Id ipgs);
    void _checkFaciesDataToGrid(Id ipgs);
    bool _keep(Id file, Id type) const;
    void _suppressAddedSamples();

  private:
    Id _numberPGS;
    std::array<Id, 2> _numberGRF{};
    std::array<Id, 2> _numberFacies{};
    Id _ngrftot;
    Id _nfactot;
    bool _flagGaus;
    bool _flagProp;
    bool _flagCheck;
    bool _flagShow;
    Id _nbtuba;
    Id _gibbsNBurn;
    Id _gibbsNIter;
    double _percent;

    RuleProp* _ruleprop; // Pointer: not to be deleted
    Model* _models[2][2]{}; // Pointer: not to be deleted

    bool _flagStat;
    bool _flagCond;
    VectorInt _iattZ{2, -1};
    Id _iptrRP;
    Id _iptrRF;
    Id _iptrDF;
    Id _iptrDN;
    Id _iptrRN;
    Id _nechInitial;
    VectorVectorBool _flagUsed{{false, false}, {false, false}};
    VectorDouble _propCst;
  };

  GSTLEARN_EXPORT Id simpgs(
    Db* dbin,
    Db* dbout,
    RuleProp* ruleprop,
    Model* model1,
    Model* model2 = nullptr,
    ANeigh* neigh = nullptr,
    Id nbsimu = 1,
    Id seed = 1321421,
    Id flag_gaus = false,
    Id flag_prop = false,
    Id flag_check = false,
    Id flag_show = false,
    Id nbtuba = 100,
    Id gibbs_nburn = 10,
    Id gibbs_niter = 100,
    double percent = 5.,
    const NamingConvention& namconv =
      NamingConvention("Facies", true, true, true, ELoc::fromKey("FACIES")));

  GSTLEARN_EXPORT Id simbipgs(
    Db* dbin,
    Db* dbout,
    RuleProp* ruleprop,
    Model* model11,
    Model* model12 = nullptr,
    Model* model21 = nullptr,
    Model* model22 = nullptr,
    ANeigh* neigh = nullptr,
    Id nbsimu = 1,
    Id seed = 43243,
    Id flag_gaus = false,
    Id flag_prop = false,
    Id flag_check = false,
    Id flag_show = false,
    Id nbtuba = 100,
    Id gibbs_nburn = 10,
    Id gibbs_niter = 100,
    double percent = 5.,
    const NamingConvention& namconv =
      NamingConvention("Facies", true, true, true, ELoc::fromKey("FACIES")));

} // namespace gstlrn
