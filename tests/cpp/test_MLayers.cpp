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
#include "API/SPDE.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/VectorHelper.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "MLayers/MLayers.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/Simulations.hpp"
#include "geoslib_f.h"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_Db-");
  Id nbtuba = 100;

  Id nech = 14;
  VectorDouble tab = {10., 10.0, 3.0, 1, 1, 10., 10.0, 4.0, 2, 1,
                      10., 10.0, 5.1, 3, 1, 10., 10.0, 6.4, 4, 1,
                      80., 10.0, 3.5, 1, 2, 81., 9.0,  3.8, 2, 2,
                      81., 9.0,  5.5, 4, 2, 12., 75.0, 2.6, 1, 3,
                      12., 75.0, 3.8, 2, 3, 12., 75.0, 4.5, 3, 3,
                      12., 75.0, 5.5, 4, 3, 65., 65.0, 3.2, 2, 4,
                      70., 70.0, 4.2, 3, 4, 80., 80.0, 6.2, 4, 4};
  auto* db = Db::createFromSamples(
    nech, ELoadBy::SAMPLE, tab, {"x1", "x2", "z1", "layer", "well"});
  db->setLocators({"x1", "x2"}, ELoc::X);
  db->setLocators({"z1"}, ELoc::Z);
  db->setLocators({"layer"}, ELoc::LAYER);
  db->display();

  auto* neigh = NeighUnique::create();
  neigh->display();

  auto* grid = DbGrid::create({100, 100});

  VectorDouble gadd = (grid->getColumn("x1") + grid->getColumn("x2")) / 200.;
  (void)grid->addColumns(gadd, "reference");

  auto* model = Model::createFromParam(ECov::SPHERICAL, 10, 0.05);
  model->setDriftIRF(0);
  (void)db->addSelectionByVariable("layer", 4, 4);
  (void)simtub(db, grid, model, neigh, 1, 0, nbtuba);
  grid->setName("Simu.z1", "bottom");
  db->clearSelection();

  auto* sills = MatrixSymmetric::createFromDiagonal({1., 3., 2., 4.});
  model =
    Model::createFromParam(ECov::CUBIC, 40., 0., 0., VectorDouble(), *sills);
  model->display();

  auto* modelT = Model::createFromParam(ECov::SPHERICAL, 10, 4);
  modelT->setMean(1000);
  (void)simtub(nullptr, grid, modelT, neigh, 1, 0, nbtuba);
  grid->setName("Simu", "Time1000");

  modelT->setMean(1200);
  (void)simtub(nullptr, grid, modelT, neigh, 1, 0, nbtuba);
  grid->setName("Simu", "Time1200");

  modelT->setMean(1400);
  (void)simtub(nullptr, grid, modelT, neigh, 1, 0, nbtuba);
  grid->setName("Simu", "Time1400");

  modelT->setMean(1600);
  (void)simtub(nullptr, grid, modelT, neigh, 1, 0, nbtuba);
  grid->setName("Simu", "Time1600");

  grid->setLocators({"Time*"}, ELoc::TIME);
  grid->display();

  Id rank = 1000;
  bool flag_same = false;
  bool flag_Z = true;
  bool flag_vel = false;
  bool flag_cumul = false;
  bool flag_ext = false;
  bool flag_std = false;
  bool flag_bayes = false;
  bool match_time = false;
  Id irf_rank = 0;

  OptDbg::setReference(rank);

  (void)multilayers_getPrior(db, grid, model, flag_same, flag_vel, flag_ext);
  (void)multilayers_kriging(
    db, grid, model, flag_same, flag_Z, flag_vel, flag_cumul, flag_ext,
    flag_bayes, flag_std, match_time, irf_rank, "reference", String(),
    "bottom");

  // MatrixDense* trace = MatrixDense::createFromVD({0, 50, 100, 0, 50, 100}, 3, 2);
  // trace->display();

  // double disc = 10.;
  // VectorDouble xp;
  // VectorDouble yp;
  // VectorDouble zp;
  // VectorDouble dd;
  // VectorDouble ddel;
  // (void)ut_trace_discretize(*trace, disc, xp, yp, dd, ddel);

  // Id checkOrder          = 2;
  // VectorVectorDouble res = interpolateVariablesToPoint(grid, xp, yp, zp, ELoc::Z,
  //                                                      checkOrder, true);
  // printVector(res, "Interpolation along lines", true, true);
  // delete trace;

  delete db;
  delete grid;
  delete neigh;
  delete model;
}
