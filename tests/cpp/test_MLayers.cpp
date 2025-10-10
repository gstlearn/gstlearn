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
#include "Basic/VectorHelper.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
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

  Id nech          = 14;
  VectorDouble tab = {10., 10., -3., 1, 1,
                      10., 10., -4., 2, 1,
                      10., 10., -5.1, 3, 1,
                      10., 10., -6.4, 4, 1,
                      80., 10., -3.5, 1, 2,
                      81., 9.5, -3.8, 2, 2,
                      81., 9.2, -5.5, 4, 2,
                      12., 75., -2.6, 1, 3,
                      12., 75., -3.8, 2, 3,
                      12., 75.1, -4.5, 3, 3,
                      12., 75.4, -5.5, 4, 3,
                      65., 65., -3.2, 2, 4,
                      70., 70., -4.2, 3, 4,
                      80., 80., -6.2, 4, 4};
  auto* db         = Db::createFromSamples(nech, ELoadBy::SAMPLE, tab,
                                           {"x1", "x2", "z1", "layer", "well"});
  db->setLocators({"x1", "x2"}, ELoc::X);
  db->setLocators({"z1"}, ELoc::Z);
  db->setLocators({"layer"}, ELoc::LAYER);
  db->display();

  auto* grid        = DbGrid::create({100, 100});
  VectorDouble gadd = VH::add(grid->getColumn("x1"), grid->getColumn("x2"));
  VH::divideConstant(gadd, 200.);
  (void)grid->addColumns(gadd, "reference");
  grid->display();

  auto* neigh = NeighUnique::create();
  neigh->display();

  auto* sills = MatrixSymmetric::createFromDiagonal({1., 3., 2., 4.});
  auto* model = Model::createFromParam(ECov::CUBIC, 40., 0., 0., VectorDouble(), *sills);
  model->display();

  // Id nbtuba    = 100;
  // auto* modelT = Model::createFromParam(ECov::SPHERICAL, 10, 4);
  // modelT->setMean(1000);
  // (void)simtub(nullptr, grid, modelT, neigh, 1, 0, nbtuba);
  // grid->setName("Simu", "Time1000");

  // modelT->setMean(1200);
  // (void)simtub(nullptr, grid, modelT, neigh, 1, 0, nbtuba);
  // grid->setName("Simu", "Time1200");

  // modelT->setMean(1400);
  // (void)simtub(nullptr, grid, modelT, neigh, 1, 0, nbtuba);
  // grid->setName("Simu", "Time1400");

  // modelT->setMean(1600);
  // (void)simtub(nullptr, grid, modelT, neigh, 1, 0, nbtuba);
  // grid->setName("Simu", "Time1600");

  // grid->setLocators({"Time*"}, ELoc::TIME);

  (void)multilayers_kriging(db, grid, model, neigh);

  MatrixDense* trace = MatrixDense::createFromVD({0, 50, 100, 0, 50, 100}, 3, 2);
  trace->display();

  double disc = 10.;
  VectorDouble xp;
  VectorDouble yp;
  VectorDouble zp;
  VectorDouble dd;
  VectorDouble ddel;
  (void)ut_trace_discretize(*trace, disc, xp, yp, dd, ddel);

  Id checkOrder          = 2;
  VectorVectorDouble res = interpolateVariablesToPoint(grid, xp, yp, zp, ELoc::Z,
                                                       checkOrder, true);
  VH::dump("Interpolation along lines", res);

  delete db;
  delete grid;
  delete neigh;
  delete model;
  delete trace;
}
