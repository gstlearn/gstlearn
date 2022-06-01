/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Model/Model.hpp"
#include "Simulation/ASimulation.hpp"
#include "Simulation/SimuRefineParam.hpp"

class Db;
class DbGrid;

class GSTLEARN_EXPORT SimuRefine: public ASimulation
{
public:
  SimuRefine(int nbsimu = 0, int seed = 4324324);
  SimuRefine(const SimuRefine &r) = delete;
  SimuRefine& operator=(const SimuRefine &r) = delete;
  virtual ~SimuRefine();

  DbGrid* simulate(DbGrid *dbin, Model *model, const SimuRefineParam& param);

private:
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

private:
  SimuRefineParam _param;
  Model* _model;
  int _ndim;
  VectorInt _nx1;
  VectorDouble _dx1;
  VectorDouble _x01;
  VectorInt _nx2;
  VectorDouble _dx2;
  VectorDouble _x02;
  int _IX[2][5];
  int _IY[2][5];
  int _IZ[2][5];
  double _XN[2][5];
  double _YN[2][5];
  double _ZN[2][5];
  double _WGT[2][2][5];
  double _STDV[2][2];
};
