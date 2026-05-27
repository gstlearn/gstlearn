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
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Boolean/ModelBoolean.hpp"
#include "Boolean/ShapeEllipsoid.hpp"
#include "Boolean/ShapeParallelepiped.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Enum/ESpaceType.hpp"
#include "Simulation/CalcSimuBoolean.hpp"
#include "Simulation/Simulations.hpp"
#include "Space/ASpaceObject.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This exercise is to demonstrate the Boolean simulation capability
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_SimBool-");

  // Global parameters
  law_set_random_seed(32131);
  Id ndim = 2;
  Id nvar = 1;
  Id nxcell = 100;
  Id nech = 100;
  bool flagConditional = true;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Generate the output grid
  VectorInt nx = {nxcell, nxcell};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Generate the data base
  auto* data = Db::createFillRandom(nech, ndim, nvar, 0, 0, 0., 0.1);
  data->getStatsAsTable().display();

  // ====================== Create Shape Dictionary ===================
  message("\n<----- Creating Shape Dictionary ----->\n");
  auto* tokens = new ModelBoolean(0.01, true);
  ShapeEllipsoid token_ellipsoid(0.4, 10., 20., 2.);
  token_ellipsoid.setFactorX2Y(1.5);
  tokens->addToken(token_ellipsoid);
  ShapeParallelepiped token_parallelepiped(0.6, 5, 7., 1.);
  tokens->addToken(token_parallelepiped);
  tokens->display();

  // ====================== Perform Boolean simulation ===================
  message("\n<----- Perform Boolean Simulation ----->\n");
  if (flagConditional)
  {
    message("- Simulation with conditioning\n");
    (void)simbool(data, grid, tokens, SimuBooleanParam(200));
  }
  else
  {
    message("- Simulation without conditioning\n");
    (void)simbool(nullptr, grid, tokens, SimuBooleanParam(200));
  }

  grid->getStatsAsTable().display();

  (void)grid->dumpToNF("grid.NF");

  delete grid;
  delete data;
  delete tokens;

  return (0);
}
