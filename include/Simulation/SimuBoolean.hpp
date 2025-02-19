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
  SimuBoolean(int nbsimu = 0, int seed = 4324324);
  SimuBoolean(const SimuBoolean &r) = delete;
  SimuBoolean& operator=(const SimuBoolean &r) = delete;
  virtual ~SimuBoolean();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int simulate(Db *dbin,
               DbGrid *dbout,
               ModelBoolean* tokens,
               const SimuBooleanParam& boolparam,
               int iptr_simu,
               int iptr_rank,
               int iptr_cover,
               bool verbose = false);

  VectorDouble extractObjects() const;

private:
  virtual bool _run() override;

  void _clearAllObjects();
  int _getNObjects(int mode = 0) const;
  int _getRankUncovered(const Db* db, int rank) const;
  int _getObjectRank(int mode, int rank);
  int _deleteObject(int mode, Db* dbin);
  static int _getAverageCount(const DbGrid* dbout,
                              const ModelBoolean* tokens,
                              const SimuBooleanParam& boolparam);
  static int _countConditioningPore(const Db* db);
  static int _countConditioningGrain(const Db* db);
  int _generatePrimary(Db* dbin,
                       DbGrid* dbout,
                       const ModelBoolean* tokens,
                       const SimuBooleanParam& boolparam,
                       bool verbose = false);
  int _generateSecondary(Db* dbin,
                         DbGrid* dbout,
                         const ModelBoolean* tokens,
                         const SimuBooleanParam& boolparam,
                         bool verbose = false);
  void _projectToGrid(DbGrid* dbout,
                      const SimuBooleanParam& boolparam,
                      int iptr_simu,
                      int iptr_rank);

private:
  std::vector<BooleanObject*> _objlist;
  mutable int _iptrCover;
};

GSTLEARN_EXPORT int simbool(Db *dbin,
                            DbGrid *dbout,
                            ModelBoolean *tokens,
                            const SimuBooleanParam& boolparam = SimuBooleanParam(),
                            int seed = 432431,
                            bool flag_simu = true,
                            bool flag_rank = true,
                            bool verbose = false,
                            const NamingConvention& namconv = NamingConvention("Boolean"));
