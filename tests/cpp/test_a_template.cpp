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

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighMoving.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  auto* mydb = Db::createFromCSV("/home/drenard//project_gstlearn/gstlearn/doc/data/Pollution/Pollution.dat");
  mydb->setLocators({"X", "Y"}, ELoc::X);
  mydb->setLocator("Zn", ELoc::Z);
  auto* mymodel = Model::createFromParam(ECov::CUBIC);
  auto* myneigh = NeighMoving::create(false, 6, 10);
  auto* mygrid  = DbGrid::createCoveringDb(mydb, VectorInt(), {0.5, 0.5}, VectorDouble(), {2, 2});
  myneigh->display();
  test_neigh(mydb, mygrid, mymodel, myneigh);

  auto* dbfmt = DbStringFormat::createFromKeys(FLAG_STATS);
  mygrid->display(dbfmt);
  return 0;
}
