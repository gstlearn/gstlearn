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

#include "Basic/AStringable.hpp"
#include "Boolean/ModelBoolean.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuBooleanParam.hpp"

namespace gstlrn
{
  class AShape;
  class BooleanObject;
  class DbGrid;
  class Db;

  /**
   * @brief Class for performing Boolean simulation
   *
   * A Boolean simulation results in drawing random object into a Field
   * conditionally to existing samples or not
   * These objects are called tokens and are generated according to:
   * - their type (extension, orientation, ...)
   * - their proportion: fix or variable
   *
   * The samples (used for conditional simulations) are defined in input Db
   * (as Z Locator variable) and are set to 0 (pore) or 1 (grain)
   * If the proportion is variable, it uses Proportion locator in output DbGrid
   */

  class GSTLEARN_EXPORT CalcSimuBoolean: public ACalcSimulation,
                                         public AStringable
  {
  public:
    CalcSimuBoolean(Id nbsimu = 0, Id seed = 4324324, bool verbose = false);
    CalcSimuBoolean(const CalcSimuBoolean& r) = delete;
    CalcSimuBoolean& operator=(const CalcSimuBoolean& r) = delete;
    virtual ~CalcSimuBoolean();

    /// Interface to AStringable
    String toString(const AStringFormat* strfmt = nullptr) const override;

    const SimuBooleanParam& getBoolparam() const { return _boolparam; }

    void setBoolParam(const SimuBooleanParam& boolparam)
    {
      _boolparam = boolparam;
    }

    void setTokens(const ModelBoolean* tokens) { _tokens = tokens; }

    void setFlagSimu(bool flag) { _flagSimu = flag; }

    void setFlagRank(bool flag) { _flagRank = flag; }

    VectorDouble extractObjects() const;

  private:
    bool _check() override;
    bool _preprocess() override;
    bool _run() override;
    bool _postprocess() override;
    void _rollback() override;

    void _clearAllObjects();
    Id _getNObjects(Id mode = 0) const;
    Id _getRankGrainUncovered(const Db* db, Id rank) const;
    Id _getObjectRank(Id mode, Id rank);
    Id _getAverageCount() const;
    Id _countConditioningPore() const;
    Id _countConditioningGrain() const;
    bool _deleteObject(Id mode, Db* dbin);
    bool _generatePrimary();
    bool _generateSecondary();
    void _projectToGrid(Id isimu);
    bool _simulate(Id isimu);
    void _resetCoverage() const;

  private:
    bool _flagSimu;
    bool _flagRank;
    SimuBooleanParam _boolparam;
    const ModelBoolean* _tokens;
    std::vector<BooleanObject> _objlist;
    mutable Id _iptrSimu;
    mutable Id _iptrRank;
    mutable Id _iptrCover;
  };

} // namespace gstlrn
