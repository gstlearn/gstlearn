#include "geoslib_f.h"

#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/ASerializable.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Model/ConsItem.hpp"
#include "Model/Constraints.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Db/Db.hpp"
#include "API/SPDE.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"

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

  model->display();

  filename = ASerializable::getTestData("Scotland","vario.ascii");
  Vario* vario = Vario::createFromNF(filename,verbose);

  vario->display();

  auto structs = {ECov::NUGGET,ECov::BESSEL_K};
  ConsItem consNug = ConsItem::define(EConsElem::SILL,0,0,0, EConsType::UPPER,0.1);
  Constraints constraints;
  constraints.addItem(&consNug);

  Option_AutoFit opt;
  OptDbg::define(EDbg::CONVERGE,1);
  int err = model->fit(vario,structs,true,opt,constraints);

  filename = ASerializable::getTestData("Scotland","model.ascii");
  model = Model::createFromNF(filename,verbose);

  NeighUnique* neighU = NeighUnique::create(ndim, false);
  neighU->display();

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Drift-");

  SPDE spde;
  if (flagSPDE)
  {
    spde.init(model,grid,temperatures,ESPDECalcMode::KRIGING);
    VectorDouble coeffs = spde.getCoeffs();
    spde.compute();
    spde.query(grid);
    grid->dumpToNF("SPDE-result.ascii",verbose);
  }
  else
  {
    kriging(temperatures, grid, model, neighU);
    grid->dumpToNF("Kriging-result.ascii",verbose);
  }

  delete temperatures;
  delete grid;
  delete model;
  delete neighU;

  return 0;
}

