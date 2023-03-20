/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Basic/OptDbg.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/File.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Model/ConsItem.hpp"
#include "Model/Constraints.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "API/SPDE.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"

int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  bool verbose = false;
  bool flagSPDE = true;
  String filename;

  filename = ASerializable::getTestData("Scotland","temperatures.ascii");
  Db* temperatures = Db::createFromNF(filename,verbose);
  temperatures->setLocator("January_temp", ELoc::Z);
  temperatures->display();

  filename = ASerializable::getTestData("Scotland","grid.ascii");
  DbGrid* grid = DbGrid::createFromNF(filename,verbose);
  grid->display();

  filename = ASerializable::getTestData("Scotland","model.ascii");
  Model* model = Model::createFromNF(filename,verbose);

  model->display();

  filename = ASerializable::getTestData("Scotland","vario.ascii");
  Vario* vario = Vario::createFromNF(filename,verbose);

  vario->display();

  auto structs = {ECov::NUGGET,ECov::BESSEL_K};
  ConsItem consNug = ConsItem::define(EConsElem::SILL,0,0,0, EConsType::UPPER,0.1);
  ConsItem consParam = ConsItem::define(EConsElem::PARAM,1, 0, 0, EConsType::EQUAL,1.);
  Constraints constraints;
  constraints.addItem(&consNug);
  constraints.addItem(&consParam);

  //OptDbg::define(EDbg::CONVERGE);
  (void) model->fit(vario,structs,constraints);
  model->display();

  NeighUnique* neighU = NeighUnique::create();
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
    (void) grid->dumpToNF("SPDE-result.ascii",verbose);
  }
  else
  {
    kriging(temperatures, grid, model, neighU);
    (void) grid->dumpToNF("Kriging-result.ascii",verbose);
  }

  DbStringFormat dbfmt(FLAG_STATS,{"*kriging"});
  grid->display(&dbfmt);

  delete temperatures;
  delete grid;
  delete model;
  delete neighU;

  return 0;
}

