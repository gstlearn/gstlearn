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
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"
#include "Enum/ESpaceType.hpp"
#include "Estimation/Estimations.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/Simulations.hpp"
#include "Simulation/SpectrumOnRN.hpp"
#include "Space/ASpaceObject.hpp"
#include "Stats/Classical.hpp"
#include "Variogram/Vario.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This test enables a comparative test of Conditional Simulation
 ** using either the Turning Band or the Spectral method.
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";

  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_SimGaussian-");

  // Global parameters
  law_set_random_seed(32131);
  Id ndim = 2;
  Id nbsimu = 3;
  Id ndat = 10;
  Id nvar = 2;
  bool flagCond = false;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Build a Model (compatible for Turning Bands and for Spectral methods)
  auto* sillExp = MatrixSymmetric::createRandomDefinitePositive(nvar);
  Model* model = Model::createFromParam(
    ECov::EXPONENTIAL, 0.2, 1., 1., VectorDouble(), *sillExp);
  auto* sillNugget = MatrixSymmetric::createRandomDefinitePositive(nvar);
  model->addCovFromParam(
    ECov::NUGGET, 0.1, 1., 1., VectorDouble(), *sillNugget);
  auto* sillGaus = MatrixSymmetric::createRandomDefinitePositive(nvar);
  model->addCovFromParam(
    ECov::GAUSSIAN, 0.15, 3., 1., VectorDouble(), *sillGaus);
  model->setMean(100.);
  model->display();

  // Build a Data set and simulate the Conditioning information
  Db* data = Db::createFillRandom(ndat, ndim, 0);
  Db* dbin = (flagCond) ? data : nullptr;

  // Simulate the conditioning information (one simulation)
  (void)simtub(
    nullptr, data, model, nullptr, 1, 5423, 100, false,
    NamingConvention("Data"));
  data->getStatsAsTable().display();

  // Generate the output grid
  Id ncell = 100;
  DbGrid* grid = DbGrid::create({ncell, ncell}, {1. / ncell, 1. / ncell});

  Id nfeatures = 100;
  // ====================== Simulation (turning bands) ====================
  message("\n<----- Simulation using Turning Bands ----->\n");
  simtub(
    dbin, grid, model, nullptr, nbsimu, 425631, nfeatures, false,
    NamingConvention("SimuTB"));
  grid->getStatsAsTable({"SimuTB*"}).display();
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    auto name = NC::getNameEncoded("SimuTB", nullptr, ivar + 1, nvar, -1);
    grid->getCorrelationAsTable({name}).display();
  }
  grid->setLocators({"SimuTB"}, ELoc::Z, 0, true);
  auto* varioTB = varioGridCalculate(grid);

  // ====================== Simulation (spectral) ====================
  message("\n<----- Simulation using Spectral Method ----->\n");
  simuSpectral(
    dbin, grid, model, nullptr, nbsimu, 425631, nfeatures, 100, false,
    NamingConvention("SimuSPT"));
  grid->getStatsAsTable({"SimuSPT*"}).display();
  for (Id ivar = 0; ivar < nvar; ivar++)
  {
    auto name = NC::getNameEncoded("SimuSPT", nullptr, ivar + 1, nvar, -1);
    grid->getCorrelationAsTable({name}).display();
  }
  grid->setLocators({"SimuSPT"}, ELoc::Z, 0, true);
  auto* varioSPT = varioGridCalculate(grid);

  // ====================== Dump into Neutral File =========================
  (void)grid->dumpToNF("Grid", EFormatNF::DEFAULT);
  (void)varioSPT->dumpToNF("Vario_SPT");
  (void)varioTB->dumpToNF("Vario_TB");

  // ====================== Free pointers ==================================
  delete data;
  delete grid;
  delete model;

  return (0);
}
