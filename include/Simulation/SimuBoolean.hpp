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
  int _getRankUncovered(const Db* db, int rank);
  int _getObjectRank(int mode, int rank);
  int _deleteObject(int mode, Db* dbin);
  int _getAverageCount(const DbGrid* dbout,
                       const ModelBoolean* tokens,
                       const SimuBooleanParam& boolparam) const;
  int _countConditioningPore(const Db* db);
  int _countConditioningGrain(const Db* db);
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
