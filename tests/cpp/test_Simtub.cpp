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
#include "Enum/EFormatNF.hpp"
#include "Enum/ESpaceType.hpp"

#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/Timer.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Space/ASpaceObject.hpp"
#include "Stats/Classical.hpp"

using namespace gstlrn;

static Db* createLocalDb(Id nech, Id ndim, Id nvar)
{
  // Coordinates
  VectorDouble tab = VH::simulateUniform(ndim * nech, 0., 50.);
  // Variable
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble tabvar = VH::simulateGaussian(nech);
    tab.insert(tab.end(), tabvar.begin(), tabvar.end());
  }

  Db* data = Db::createFromSamples(nech, ELoadBy::COLUMN, tab);
  data->setNameByUID(1, "x1");
  data->setNameByUID(2, "x2");

  data->setLocatorByUID(1, ELoc::X, 0);
  data->setLocatorByUID(2, ELoc::X, 1);

  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    data->setNameByUID(3 + ivar, "Var");
    data->setLocatorByUID(3 + ivar, ELoc::Z, ivar);
  }
  return data;
}

/**
 * Test of heterotopic kriging and simulation
 * Returns 1 if the rest of the test must be discarded
 */
Id st_mini_test()
{
  // Global parameters
  Id nbsimu     = 2;
  Id mode       = 0; // 1: NCsimu; 2: Kriging; 3: CDsimu; 0: All
  bool debug    = true;
  bool end_test = false;
  OptCst::define(ECst::NTCOL, -1);

  // Create the Data Base
  Db* db = Db::createFillRandom(5, 2, 2);
  db->setValue("z-1", 0, TEST);
  db->setValue("z-1", 2, TEST);
  db->setValue("z-1", 4, TEST);

  db->setValue("z-2", 1, TEST);
  db->setValue("z-2", 2, TEST);

  VectorDouble means = {100., 0};

  // Modify the variables by adding their mean
  VectorDouble z1 = db->getColumn("z-1");
  VH::addConstant(z1, means[0]);
  db->setColumn(z1, "z-1");

  VectorDouble z2 = db->getColumn("z-2");
  VH::addConstant(z2, means[1]);
  db->setColumn(z2, "z-2");

  DbStringFormat* dbfmt = DbStringFormat::createFromFlags(false, false, false, false, true);
  db->display(dbfmt);

  MatrixSymmetric* sills = MatrixSymmetric::createFromVD({3, 1, 1, 2});
  Model* model           = Model::createFromParam(ECov::SPHERICAL, 1, 1, 1, VectorDouble(), *sills);
  delete sills;
  model->setMeans(means);

  NeighMoving* neigh = NeighMoving::create(false, 100, 10);

  // Creating the output grid
  DbGrid* grid = DbGrid::create({2, 2});

  // Set the debugging option
  if (debug) OptDbg::setReference(1);

  // Perform non-conditional simulations
  if (mode == 0 || mode == 1)
    (void)simtub(nullptr, grid, model, neigh, nbsimu);

  // Perform Kriging
  if (mode == 0 || mode == 2)
    (void)kriging(db, grid, model, neigh);

  // Perform conditional simulations
  if (mode == 0 || mode == 3)
    (void)simtub(db, grid, model, neigh, nbsimu);

  grid->display(dbfmt);

  OptDbg::setReference(0);

  delete db;
  delete model;
  delete grid;

  return end_test;
}

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";

  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_Simtub-");

  // Global parameters
  law_set_random_seed(32131);
  bool verbose     = 1;
  Id ndim          = 2;
  Id nvar          = 1;
  Id nbsimu        = 3;
  DbGrid* grid_res = nullptr;
  defineDefaultSpace(ESpaceType::RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS, {"Simu*"});

  // Perform a preliminary test to check heterotopic conditional simulation
  if (st_mini_test()) return 0;

  // Generate the output grid
  VectorInt nx = {50, 50};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Generate the data base
  Id nech  = 100;
  Db* data = createLocalDb(nech, ndim, nvar);
  data->display(&dbfmt);

  // Create the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 10., 1.);
  model->display();

  // Creating a Moving Neighborhood
  NeighMoving* neighM = NeighMoving::create(false, 25);
  neighM->display();

  // Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create();
  neighU->display();

  Timer timer;
  // ====================== Simulation (turning bands) ====================
  message("\n<----- Simulation (Moving Neighborhood) ----->\n");
  delete grid_res;
  grid_res = grid->clone();
  simtub(data, grid_res, model, neighM, nbsimu);
  grid_res->display(&dbfmt);
  (void)grid_res->dumpToNF("Moving.NF", EFormatNF::DEFAULT, verbose);

  message("\n<----- Simulation (Unique Neighborhood) ----->\n");
  delete grid_res;
  grid_res = grid->clone();
  simtub(data, grid_res, model, neighU, nbsimu);
  grid_res->display(&dbfmt);
  (void)grid_res->dumpToNF("Unique.NF", EFormatNF::DEFAULT, verbose);

  timer.displayIntervalMilliseconds("Turning Band Simulations", 773);

  // ====================== Free pointers ==================================
  delete neighM;
  delete neighU;
  delete data;
  delete grid;
  delete grid_res;
  delete model;

  return (0);
}
