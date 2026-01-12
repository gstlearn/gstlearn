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
#include "Basic/OptCustom.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CorGaussianMixture.hpp"
#include "Covariances/CorGneiting.hpp"
#include "Covariances/CovContext.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"

using namespace gstlrn;
/**
 * This file is meant to test Kriging with Gneiting Model
 */
int main(int argc, char* argv[])
{
  OptCustom::define("ompthreads", 1);
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  auto space1d = SpaceRN::create(1);
  auto space2d = SpaceRN::create(2);
  auto sp      = SpaceComposite::create({space2d, space1d});
  sp->display();
  setDefaultSpace(sp);

  // creating the time trace covariance
  CovContext ctxt1d(1, space1d);
  double alpha  = 1.0;
  double beta   = 1.0;
  double scaleT = 5.3;
  auto* corT    = new CorAniso(CovContext(1, 1), ECov::CAUCHY_GEN);
  corT->setParam(alpha, 0);        // alpha in (0,2]
  corT->setParam(beta * 2 / 2, 1); // beta*d/2 with beta in (0,1]
  corT->setScaleDim(0, scaleT);

  // creating the space trace covariance
  CovContext ctxt2d(1, space2d);
  VectorDouble scales = {2., 3.};
  VectorDouble params = {0.5};
  VectorDouble kappas = {1.0};
  VectorDouble angles = {30.0, 0.0};
  auto* corS          = CorGaussianMixture::create(
    ctxt2d,
    ECov::MATERN,
    params,
    kappas,
    scales,
    angles,
    false);
  double sep = 1.;
  CorGneiting corGneiting(corS, corT, sep);
  message("Space dimension of Gneiting Covariance = %d\n", corGneiting.getNDim());

  // Testing the covariance calculation between two points
  VectorDouble coords1 = {12., 3., 1.};
  VectorDouble coords2 = {4., 5., 2.};
  SpacePoint p1(sp);
  SpacePoint p2(sp);
  p1.setCoords(coords1);
  p2.setCoords(coords2);
  double cres = corGneiting.evalCov(p1, p2);
  std::cout << "Value of Gneiting (by Covariance) = " << cres << std::endl;

  // Create the Data Base
  Id ndim  = 3;
  Id ndat  = 10;
  Id nvar  = 1;
  Db* data = Db::createFillRandom(ndat, ndim, nvar);

  // Create the Target
  VectorInt nx    = {5, 5, 2};
  VectorDouble dx = {1. / nx[0], 1. / nx[1], 1. / nx[2]};
  DbGrid* grid    = DbGrid::create(nx, dx);

  // Create the Model
  auto* model = new ModelGeneric();
  model->setCov(&corGneiting);
  model->evalCov(p1, p2);
  message("Model dimension = %d\n", model->getNDim());
  std::cout << "Value of Gneiting (by Model) = " << cres << std::endl;

  // Create the Unique neighborhood
  NeighUnique* neigh = NeighUnique::create(false, sp);
  // Launch Kriging
  (void)kriging(data, grid, model, neigh);

  // Display a summary of the results
  DbStringFormat dbfmtKriging(FLAG_STATS);
  grid->display(&dbfmtKriging);

  delete corT;
  delete corS;
  delete data;
  delete grid;
  delete neigh;
  delete model;

  return (0);
}
