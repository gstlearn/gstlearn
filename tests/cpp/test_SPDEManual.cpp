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
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "LinearOp/InvNuggetOp.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/SPDEOp.hpp"
#include "LinearOp/ShiftOpMatrix.hpp"
#include "Mesh/AMesh.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"

#define __USE_MATH_DEFINES
#include <cmath>

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_SPDEManual-");
  Id seed = 10355;
  bool verbose = true;
  law_set_random_seed(seed);

  ///////////////////////
  // Creating the Db Grid
  VectorInt nx = {101, 101};
  DbGrid* grid = DbGrid::create(nx);
  if (verbose) grid->display();

  //////////////////////
  // Creating the Mesh
  MeshETurbo mesh(grid);
  if (verbose) mesh.display();

  ///////////////////////
  // Creating the Model
  auto* model = Model::createFromParam(ECov::MATERN, 1., 1., 1., {10., 45.});
  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  model->getCovAniso(0)->makeAngleNoStatFunctional(&spirale);
  if (verbose) model->display();

  //////////////////////////////////
  // Creating the Precision Operator
  auto* cova = model->getCovAniso(0);
  ShiftOpMatrix S(&mesh, cova, grid);
  PrecisionOp Qsimu(&S, cova, verbose);

  /////////////////////////
  // Simulation (Chebyshev)
  VectorDouble gridSimu = Qsimu.simulateOne();
  grid->addColumns(gridSimu, "SimuNC", ELoc::Z);

  ///////////////////////////
  // Creating Data
  auto ndata = 1000;
  auto* dat = Db::createFromBox(
    ndata, grid->getCoorMinimum(), grid->getCoorMaximum(), 432432);

  ///////////////////////////////////////////////
  // Non-conditional Simulation at Data locations
  ProjMatrix B(dat, &mesh);
  VectorDouble datval(ndata);
  B.mesh2point(gridSimu, datval);
  dat->addColumns(datval, "SimuNC", ELoc::Z);

  //////////
  // Kriging
  auto napices = S.getSize();
  VectorDouble rhs(napices);
  B.point2mesh(datval, rhs);

  double nug = 0.1;
  for (Id i = 0; i < napices; i++) rhs[i] /= nug;

  double vardata = 0.01;
  PrecisionOp Qkriging(&S, cova);

  auto* modelNug = Model::createFromParam(ECov::NUGGET, 0., vardata);
  InvNuggetOp Bnug(dat, modelNug);
  SPDEOp A(&Qkriging, &B, &Bnug);

  VectorDouble vc(napices);

  vc = A.evalInverse(rhs);
  grid->addColumns(vc, "Kriging");

  // Statistics
  DbStringFormat dsf(FLAG_RESUME | FLAG_STATS);
  dat->display(&dsf);
  grid->display(&dsf);
  (void)grid->dumpToNF("spde.NF");

  delete dat;
  delete grid;
  delete model;
  delete modelNug;
  return 0;
}
