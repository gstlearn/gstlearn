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
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CorGneiting.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpaceComposite.hpp"
#include "Estimation/CalcKriging.hpp"
/**
 * This file is meant to test Kriging with Gneiting Model
 */
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  auto space1d = SpaceRN::create(1);
  auto space2d = SpaceRN::create(2);
  auto sp      = SpaceComposite::create({space2d, space1d});
  sp->display();
  setDefaultSpace(sp);

  double scaleT  = 5.3;
  CovAniso* covT = CovAniso::createFromParam(ECov::EXPONENTIAL,
                                             scaleT,
                                             1.,
                                             1.,
                                             VectorDouble(),
                                             MatrixSquareSymmetric(),
                                             VectorDouble(),
                                             space1d,
                                             false);

  VectorDouble scales = {2., 3.};
  CovAniso* covS      = CovAniso::createFromParam(ECov::EXPONENTIAL,
                                                  1.,
                                                  1.,
                                                  1.,
                                                  scales,
                                                  MatrixSquareSymmetric(),
                                                  VectorDouble(),
                                                  space2d,
                                                  false);

  double sep           = 1.;
  CorGneiting covGneiting = CorGneiting(covS->getCorAniso(), covT->getCorAniso(), sep);
  message("Space dimension of Gneiting Covariance = %d\n", covGneiting.getNDim());

  // Testing the covariance calculation between two points
  VectorDouble coords1 = {12., 3., 1.};
  VectorDouble coords2 = { 4., 5., 2.};
  SpacePoint p1(sp);
  SpacePoint p2(sp);
  p1.setCoords(coords1);
  p2.setCoords(coords2);
  double cres = covGneiting.evalCov(p1,p2);
  std::cout << "Value of Gneiting (by Covariance) = " << cres <<std::endl;

  // Create the Data Base
  int ndim = 3;
  int ndat = 10;
  int nvar = 1;
  Db* data = Db::createFillRandom(ndat, ndim, nvar);

  // Create the Target
  VectorInt nx    = {5, 5, 2};
  VectorDouble dx = {1. / nx[0], 1. / nx[1], 1. / nx[2]};
  DbGrid* grid = DbGrid::create(nx, dx);

  // Create the Model
  ModelGeneric* model = new ModelGeneric();
  model->setCov(&covGneiting);
  model->evalCov(p1, p2);
  message("Model dimension = %d\n", model->getNDim());
  std::cout << "Value of Gneiting (by Model) = " << cres << std::endl;

  // Create the Unique neighborhood
  NeighUnique* neigh = NeighUnique::create(false, sp);
  // Launch Kriging
  (void) kriging(data, grid, model, neigh);

  // Display a summary of the results
  DbStringFormat dbfmtKriging(FLAG_STATS);
  grid->display(&dbfmtKriging);

  delete covT;
  delete covS;
  delete data;
  delete grid;
  delete neigh;
  delete model;

  return(0);
}
