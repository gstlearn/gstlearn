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
#include "Model/Model.hpp"
#include "MultiLayers/MGibbs.hpp"
#include "Polygon/Polygons.hpp"
#include "geoslib_define.h"

using namespace gstlrn;
/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the MultiLayers 2.5D estimation / simulation
 ** with constraints (MGibbs)
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-");

  Id nlayer = 5;
  Id nline  = 10;
  auto* db  = Db::createVerticalWellsFillRandom(nline, nlayer, 0., 1. / nlayer, 0.25);

  auto* dbfmt = DbStringFormat::createFromFlags(false, false, false, false, true);
  db->display(dbfmt);

  auto* dbgrid = DbGrid::create({100, 100}, {0.01, 0.01});
  dbgrid->display();

  auto* model = Model::createFromParam(ECov::MATERN, 0.3, 1., 2.);
  model->display();

  Id niter     = 10;
  Id nbsimu    = 1;
  Id seed      = 13243;
  bool verbose = true;
  (void)MLayers_spde(db, dbgrid, model, nlayer, niter, nbsimu, -1,
                     false, false, false, false, seed, verbose);

  dbgrid->display();
  dbgrid->dumpToNF("MGrid2D");

  delete dbgrid;
  delete db;
  delete dbfmt;
  return 0;
}
