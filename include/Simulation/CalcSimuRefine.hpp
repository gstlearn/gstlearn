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

class Db;
class DbGrid;

class GSTLEARN_EXPORT CalcSimuRefine: public ACalcSimulation
{
public:
  CalcSimuRefine(int nbsimu = 0, int seed = 4324324);
  CalcSimuRefine(const CalcSimuRefine &r) = delete;
  CalcSimuRefine& operator=(const CalcSimuRefine &r) = delete;
  virtual ~CalcSimuRefine();

  const SimuRefineParam& getParam() const { return _param; }
  void setParam(const SimuRefineParam& param) { _param = param; }
  DbGrid* getResultingGrid() const { return _dbres->clone(); }

private:
  virtual bool _check() override;
  virtual bool _preprocess() override;
  virtual bool _run() override;

  void _dim_1_to_2(DbGrid *db);
  void _dim_2_to_1(DbGrid *db);
  int _kriging_define();
  void _neigh_simfine(int type, int rank, int idx, int idy, int idz);
  void _merge_data(DbGrid *db1, int iatt1, DbGrid *db2, int iatt2);
  double _read(DbGrid* db, int iatt, int ix0, int iy0, int iz0, int idx, int idy, int idz);
  static void _write(DbGrid* db, int iatt, int ix0, int iy0, int iz0, double value);
  void _truncate_result(DbGrid *db2, int iatt2, DbGrid *db1, int iatt1);
  int _kriging_solve(int type,
                     int rank,
                     int nb,
                     bool verbose = false);
  void _simulate_nodes(DbGrid *db, int iatt);
  void _simulate_target(DbGrid* db, int type, int iatt, int ix0, int iy0, int iz0);
  int _simulate();

private:
  SimuRefineParam _param;
  VectorInt _nx1;
  VectorDouble _dx1;
  VectorDouble _x01;
  VectorInt _nx2;
  VectorDouble _dx2;
  VectorDouble _x02;

  DbGrid* _dbres; // Resulting Grid
  int _IXYZ[3][2][5];
  double _XYZN[3][2][5];
  double _WGT[2][2][5];
  double _STDV[2][2];
};

GSTLEARN_EXPORT DbGrid* simulation_refine(DbGrid* dbin,
                                          Model* model,
                                          const SimuRefineParam& param,
                                          int seed                        = 432432,
                                          const NamingConvention& namconv = NamingConvention("Refine"));
