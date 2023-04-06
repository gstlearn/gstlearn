/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "geoslib_define.h"

#include "Model/Model.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuRefineParam.hpp"

class Db;
class DbGrid;

class GSTLEARN_EXPORT SimuRefine: public ACalcSimulation
{
public:
  SimuRefine(int nbsimu = 0, int seed = 4324324);
  SimuRefine(const SimuRefine &r) = delete;
  SimuRefine& operator=(const SimuRefine &r) = delete;
  virtual ~SimuRefine();

  DbGrid* simulate(DbGrid *dbin, Model *model, const SimuRefineParam& param);

private:
  virtual bool _run() override;

  void _dim_1_to_2(DbGrid *db);
  void _dim_2_to_1(DbGrid *db);
  int _kriging_define();
  void _neigh_simfine(int type, int rank, int idx, int idy, int idz);
  void _merge_data(DbGrid *db1, int iatt1, DbGrid *db2, int iatt2);
  double _read(DbGrid *db,
               int iatt,
               int ix0,
               int iy0,
               int iz0,
               int idx,
               int idy,
               int idz);
  void _write(DbGrid *db, int iatt, int ix0, int iy0, int iz0, double value);
  void _truncate_result(DbGrid *db2, int iatt2, DbGrid *db1, int iatt1);
  int _kriging_solve(int type,
                     int rank,
                     int nb,
                     bool verbose = false);
  void _simulate_nodes(DbGrid *db, int iatt);
  void _simulate_target(DbGrid *db,
                        int type,
                        int iatt,
                        int ix0,
                        int iy0,
                        int iz0);
  int _getNDim() const;

private:
  SimuRefineParam _param;
  Model* _model;
  VectorInt _nx1;
  VectorDouble _dx1;
  VectorDouble _x01;
  VectorInt _nx2;
  VectorDouble _dx2;
  VectorDouble _x02;
  int _IXYZ[3][2][5];
  double _XYZN[3][2][5];
  double _WGT[2][2][5];
  double _STDV[2][2];
};
