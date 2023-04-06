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
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include "Enum/ELoc.hpp"

#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/File.hpp"
#include "Basic/Indirection.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Indirection class
 ** This class is instantiated as a member of MeshETurbo class when a
 ** selection is present in the input Grid
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  DbStringFormat dbfmt(FLAG_STATS);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Indirect-");

  // Creating the small Grid Db
  // All characteristics along X and Y are different (on purpose)
  DbGrid* grid = DbGrid::create({5,4}, {1.,2.}, {10.,20.});
  int nech = grid->getSampleNumber();

  // Add a selection to the Grid
  VectorDouble x1 = grid->getColumn("x1");
  VectorDouble sel;
  for (int i = 0; i < nech; i++)
  {
    double value = (x1[i] > 12.) ? 1. : 0.;
    sel.push_back(value);
  }
  grid->addColumns(sel,"sel",ELoc::SEL);
  grid->display();

  // Creating a Mesh based on the Masked Grid (in order to start Indirections)

  MeshETurbo mesh = MeshETurbo(grid, false, false, 1);
  mesh.display();
  mesh.printMeshes(2);

  // Check back-and-forth between relative and absolute
  Indirection indirect = mesh.getGridIndirect();
  int igrid_rel = 3;
  message("Starting from the Relative grid node = %d\n",igrid_rel);
  int igrid_abs = indirect.getRToA(igrid_rel);
  message("Corresponding absolute grid node = %d\n",igrid_abs);
  int igrid_rel2 = indirect.getAToR(igrid_abs);
  message("Ending Relative grid node = %d\n",igrid_rel2);

  // Performing the same task with the Indirection stored in Integer Arrays
  // This should be quicker but more space consuming

  MeshETurbo mesh0 = MeshETurbo(grid, false, false, 0);
  mesh0.display();
  mesh0.printMeshes(2);

  delete grid;
  return 0;
}

