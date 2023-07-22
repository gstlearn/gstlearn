/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_d.h"

#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Simulation/SimuSubstitutionParam.hpp"
#include "Simulation/CalcSimuSubstitution.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This exercise is to demonstrate the Substitution simulation capability
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("SimSub-");

  // Global parameters
  int ndim = 2;
  int seed = 3322;
  int nxcell = 100;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Generate the output grid
  VectorInt nx = {nxcell,nxcell};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // ====================== Create Parameter File ===================
  message("\n<----- Creating Parameter File ----->\n");
  int nfacies = 3;
  double intensity = 1.;
  SimuSubstitutionParam subparam(nfacies, intensity);
  subparam.display();

  // ====================== Perform Boolean simulation ===================
  message("\n<----- Perform Substitution Simulation ----->\n");
  (void) substitution(grid, subparam, seed, false);

  subparam.setFactor(0.6);
  subparam.setVector({0.,1.});
  subparam.setFlagOrient(true);
  (void) substitution(grid, subparam, seed, false);

  (void) grid->dumpToNF("grid.ascii");

  delete grid;

  return (0);
}
