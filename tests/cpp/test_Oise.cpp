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
#include "geoslib_define.h"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbHelper.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"
#include "utils.hpp"
#include "Basic/Timer.hpp"

using namespace gstlrn;

/**
 * This test is meant to execute Oise Kriging in CPP (heaptrack diagnostic)
 */


Db* createOiseData(bool selvario = false)
{
  message("Creating data...\n");
  String filename = getTestData("Alluvial","Oise_Data_ThicknessSides.ascii");
  Db* data = Db::createFromNF(filename);
  data->clearLocators(ELoc::L);
  data->clearLocators(ELoc::U);
  VectorDouble error(data->getNSample(), 0.1);
  data->addColumns(error, "error");
  data->setLocator("error", ELoc::V);
  data->clearLocators(ELoc::SEL);
  if (selvario)
  {
    VectorDouble thicknessSides = data->getColumn("ThicknessSides");
    VectorDouble sel;
    for (size_t idat = 0; idat < thicknessSides.size(); ++idat)
    {
      if (thicknessSides[idat] != 0.)
        sel.push_back(1.);
      else
        sel.push_back(0.);
    }
    data->addColumns(sel, "selVario");
    data->setLocator("selVario", ELoc::SEL);
  }
  data->setLocator("ThicknessSides", ELoc::Z, 0, true);
  data->display();
  return data;
}

DbGrid* createOiseGrid()
{
  message("\nCreating grid...\n");
  String filename = getTestData("Alluvial", "Oise_GridAnglesModifFinal.ascii");
  DbGrid* grid = DbGrid::createFromNF(filename);
  if (!grid)
  {
    messerr("Cannot create grid from file %s", filename.c_str());
    messerr("Please add a simlink to the required file in doc/data/Alluvial/");
    return nullptr;
  }
  grid->setLocator("Polygon.*", ELoc::SEL);
  grid->deleteColumns({"spde*"});
  VectorDouble ranges(grid->getNSample(), 200.);
  grid->addColumns(ranges, "ranges0");
  grid->addColumns(ranges, "ranges1");
  grid->display();
  return grid;
}

Model* makeOiseModel(const DbGrid* grid = nullptr, bool stationary = false, Id drift = 1)
{
  message("\nCreating model...\n");
  VectorDouble ranges = {8000,800};
  VectorDouble angles = {10,0};
  Model* model = Model::createFromParam(ECov::MATERN, 1., 1., 1., ranges, MatrixSymmetric(), angles);
  if (grid != nullptr && !stationary)
      model->getCovAniso(0)->makeAngleNoStatDb("Migrate.angles_interp", 0, grid);
  else
      model->getCovAniso(0)->setAnisoAngle(0,45);
  model->setDriftIRF(drift);
  model->display();
  return model;
}

int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  bool stationary = false;
  Db* data = createOiseData(true);
  DbGrid* grid = createOiseGrid();
  Model* model = makeOiseModel(grid, stationary, 1);

  // Kriging (stationary)
  message("\nKriging (non-stationary)...\n");
  Timer timer;
  NeighUnique* neigh = NeighUnique::create();
  Id error = kriging(data, grid, model, neigh);
  message("Kriging (non-stationary) done - status=%d.\n", error);
  timer.displayIntervalSeconds("Kriging (non-stationary)", 5310);

  delete data;
  delete grid;
  delete model;
  delete neigh;

  return 0;
}
