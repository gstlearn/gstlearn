#include "geoslib_f.h"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Neigh/Neigh.hpp"
#include "Basic/ASerializable.hpp"

int main(int /*argc*/, char */*argv*/[])
{
  bool verbose = false;
  bool flagSPDE = true;
  int ndim = 2;
  String filename;

  filename = ASerializable::getTestData("Scotland","temperatures.ascii");
  Db* temperatures = Db::createFromNF(filename,false,verbose);
  temperatures->setLocator("January_temp", ELoc::Z);
  temperatures->display();

  filename = ASerializable::getTestData("Scotland","grid.ascii");
  Db* grid = Db::createFromNF(filename,true,verbose);
  grid->display();

  filename = ASerializable::getTestData("Scotland","model.ascii");
  Model* model = Model::createFromNF(filename,verbose);

//  CovContext ctxt(1,ndim);
//  Model* model = Model::create(ctxt);
//  CovLMC covs(ctxt.getSpace());
//  double range = 161.;
//  double param = 1.;
//  double sill = 1.119260;
//  CovAniso cova = CovAniso(ECov::BESSEL_K,range,param,sill,ctxt);
//  covs.addCov(&cova);
//  double sill_nug = 0.131347;
//  cova = CovAniso(ECov::NUGGET,0.,0.,sill_nug,ctxt);
//  covs.addCov(&cova);
//
//  model->setCovList(&covs);

  model->display();

  Neigh* neigh = Neigh::createUnique(ndim);
  neigh->display();

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Drift-");

  SPDE spde;
  if (flagSPDE)
  {
    spde.init(model,grid,temperatures,ESPDECalcMode::KRIGING);
    spde.compute();
    spde.query(grid);
    grid->dumpToNF("SPDE-result.ascii",verbose);
  }
  else
  {
    kriging(temperatures, grid, model, neigh);
    grid->dumpToNF("Kriging-result.ascii",verbose);
  }

  delete temperatures;
  delete grid;
  delete model;

  return 0;
}

