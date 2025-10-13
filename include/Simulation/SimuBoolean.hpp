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

#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuBooleanParam.hpp"
#include "Boolean/ModelBoolean.hpp"
#include "Basic/AStringable.hpp"

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


class GSTLEARN_EXPORT SimuBoolean: public ACalcSimulation, public AStringable
{
public:
  SimuBoolean(Id nbsimu = 0, Id seed = 4324324);
  SimuBoolean(const SimuBoolean &r) = delete;
  SimuBoolean& operator=(const SimuBoolean &r) = delete;
  virtual ~SimuBoolean();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  Id simulate(Db *dbin,
               DbGrid *dbout,
               ModelBoolean* tokens,
               const SimuBooleanParam& boolparam,
               Id iptr_simu,
               Id iptr_rank,
               Id iptr_cover,
               bool verbose = false);

  VectorDouble extractObjects() const;

private:
  bool _run() override;

  void _clearAllObjects();
  Id _getNObjects(Id mode = 0) const;
  Id _getRankUncovered(const Db* db, Id rank) const;
  Id _getObjectRank(Id mode, Id rank);
  Id _deleteObject(Id mode, Db* dbin);
  static Id _getAverageCount(const DbGrid* dbout,
                              const ModelBoolean* tokens,
                              const SimuBooleanParam& boolparam);
  static Id _countConditioningPore(const Db* db);
  static Id _countConditioningGrain(const Db* db);
  Id _generatePrimary(Db* dbin,
                       DbGrid* dbout,
                       const ModelBoolean* tokens,
                       const SimuBooleanParam& boolparam,
                       bool verbose = false);
  Id _generateSecondary(Db* dbin,
                         DbGrid* dbout,
                         const ModelBoolean* tokens,
                         const SimuBooleanParam& boolparam,
                         bool verbose = false);
  void _projectToGrid(DbGrid* dbout,
                      const SimuBooleanParam& boolparam,
                      Id iptr_simu,
                      Id iptr_rank);

private:
  std::vector<BooleanObject*> _objlist;
  mutable Id _iptrCover;
};

GSTLEARN_EXPORT Id simbool(Db *dbin,
                            DbGrid *dbout,
                            ModelBoolean *tokens,
                            const SimuBooleanParam& boolparam = SimuBooleanParam(),
                            Id seed = 432431,
                            bool flag_simu = true,
                            bool flag_rank = true,
                            bool verbose = false,
                            const NamingConvention& namconv = NamingConvention("Boolean"));
}
