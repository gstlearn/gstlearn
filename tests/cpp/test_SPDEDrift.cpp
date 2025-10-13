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
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Enum/EFormatNF.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Model/ConsItem.hpp"
#include "Model/Constraints.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Variogram/Vario.hpp"
#include "utils.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  bool verbose = true;
  // This crashes under MingGW/Windows due to compatibility issue with getTestData
  // ASerializable::setPrefixName("test_SPDEDrift-");

  String filename;

  filename         = getTestData("Scotland", "temperatures.ascii");
  Db* temperatures = Db::createFromNF(filename, verbose);
  temperatures->setLocator("January_temp", ELoc::Z, 0);
  temperatures->display();

  filename     = getTestData("Scotland", "grid.ascii");
  DbGrid* grid = DbGrid::createFromNF(filename, verbose);
  grid->display();

  filename     = getTestData("Scotland", "model.ascii");
  Model* model = Model::createFromNF(filename, verbose);

  model->display();

  filename     = getTestData("Scotland", "vario.ascii");
  Vario* vario = Vario::createFromNF(filename, verbose);

  vario->display();

  auto structs       = {ECov::NUGGET, ECov::MATERN};
  ConsItem consNug   = ConsItem::define(EConsElem::SILL, 0, 0, 0, EConsType::UPPER, 0.1);
  ConsItem consParam = ConsItem::define(EConsElem::PARAM, 1, 0, 0, EConsType::EQUAL, 1.);
  Constraints constraints;
  constraints.addItem(&consNug);
  constraints.addItem(&consParam);

  // OptDbg::define(EDbg::CONVERGE);
  (void)model->fit(vario, structs, constraints);
  model->display();

  NeighUnique* neighU = NeighUnique::create();
  neighU->display();

  //////////////////////
  /// Kriging using SPDE
  Id useCholesky = 0;
  krigingSPDE(temperatures, grid, model, true, false, useCholesky);

  ///////////////////////
  /// Traditional Kriging
  kriging(temperatures, grid, model, neighU);

  (void)grid->dumpToNF("Grid.ascii", EFormatNF::DEFAULT, verbose);

  DbStringFormat dbfmt(FLAG_STATS, {"*estim"});
  grid->display(&dbfmt);

  delete temperatures;
  delete grid;
  delete model;
  delete neighU;
  delete vario;

  return 0;
}
