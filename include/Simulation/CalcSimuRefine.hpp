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

#include "Model/Model.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuRefineParam.hpp"

namespace gstlrn
{
class Db;
class DbGrid;
class GSTLEARN_EXPORT CalcSimuRefine: public ACalcSimulation
{
public:
  CalcSimuRefine(Id nbsimu = 0, Id seed = 4324324);
  CalcSimuRefine(const CalcSimuRefine& r)            = delete;
  CalcSimuRefine& operator=(const CalcSimuRefine& r) = delete;
  virtual ~CalcSimuRefine();

  const SimuRefineParam& getParam() const { return _param; }
  void setParam(const SimuRefineParam& param) { _param = param; }
  DbGrid* getResultingGrid() const { return _dbres->clone(); }

private:
  bool _check() override;
  bool _preprocess() override;
  bool _run() override;

  void _dim_1_to_2(DbGrid* db);
  void _dim_2_to_1(DbGrid* db);
  Id _kriging_define();
  void _neigh_simfine(Id type, Id rank, Id idx, Id idy, Id idz);
  void _merge_data(DbGrid* db1, Id iatt1, DbGrid* db2, Id iatt2);
  double _read(DbGrid* db, Id iatt, Id ix0, Id iy0, Id iz0, Id idx, Id idy, Id idz);
  static void _write(DbGrid* db, Id iatt, Id ix0, Id iy0, Id iz0, double value);
  void _truncate_result(DbGrid* db2, Id iatt2, DbGrid* db1, Id iatt1);
  Id _kriging_solve(Id type,
                     Id rank,
                     Id nb,
                     bool verbose = false);
  void _simulate_nodes(DbGrid* db, Id iatt);
  void _simulate_target(DbGrid* db, Id type, Id iatt, Id ix0, Id iy0, Id iz0);
  Id _simulate();

private:
  SimuRefineParam _param;
  VectorInt _nx1;
  VectorDouble _dx1;
  VectorDouble _x01;
  VectorInt _nx2;
  VectorDouble _dx2;
  VectorDouble _x02;

  DbGrid* _dbres; // Resulting Grid
  Id _IXYZ[3][2][5];
  double _XYZN[3][2][5];
  double _WGT[2][2][5];
  double _STDV[2][2];
};

GSTLEARN_EXPORT DbGrid* simulation_refine(DbGrid* dbin,
                                          Model* model,
                                          const SimuRefineParam& param,
                                          Id seed                        = 432432,
                                          const NamingConvention& namconv = NamingConvention("Refine"));
} // namespace gstlrn