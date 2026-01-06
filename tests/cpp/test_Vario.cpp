/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/*                                                                            */
/* This file is meant to demonstrate the process of using PGS                 */
/*                                                                            */
/******************************************************************************/
#include "Basic/File.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Enum/ECalcVario.hpp"
#include "Enum/ECov.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Variogram/VMap.hpp"
#include "Variogram/Vario.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  Id error = 1;
  Id ndim  = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  CovContext ctxt(1, 2, 1.); // use default space

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  Id nech = 1000;
  Db* db  = Db::createFromBox(nech, {0., 0.}, {1., 1.}, 3242);
  db->display();

  // Creating a grid covering the same space
  VectorInt nx    = {100, 100};
  VectorDouble dx = {0.01, 0.01};
  DbGrid* grid    = DbGrid::create(nx, dx);
  grid->display();

  // Creating the Model(s) of the Underlying GRF(s)
  Model models(ctxt);
  CovAnisoList covs(ctxt);
  double range1 = 0.2;
  CovAniso cova1(ECov::MATERN, range1, 1., 1., ctxt);
  covs.addCov(cova1);
  models.setCovAnisoList(&covs);
  models.display();

  // Perform a non-conditional simulation on the Db and on the Grid
  error = simtub(nullptr, db, &models);
  db->display();
  error = simtub(nullptr, grid, &models);
  grid->display();

  // ===============
  // On Data samples
  // ===============

  mestitle(0, "Experimental variogram on Data Samples");
  Id nlag                 = 20;
  VarioParam* varioparamP = VarioParam::createMultiple(2, nlag, 0.5 / nlag);
  Vario* variop           = Vario::computeFromDb(*varioparamP, db, ECalcVario::VARIOGRAM);
  variop->display();
  delete varioparamP;

  // Fitting the experimental variogram of Underlying GRF (with constraint that total sill is 1)
  Model model(ctxt);
  VectorECov covas {ECov::MATERN, ECov::EXPONENTIAL};
  model.fit(variop, covas, false);
  model.display();
  delete variop;

  // With a quick entry
  mestitle(0, "Experimental variogram on Data Samples (quick entry)");
  auto* variop2 = variogramCalculate(db, ECalcVario::VARIOGRAM, true, 5, 0.1);
  variop2->display();
  delete variop2;

  // ===============
  // On Grid samples
  // ===============

  mestitle(0, "Experimental variogram on Grid");
  VarioParam* varioparamG = VarioParam::createMultipleFromGrid(grid, nlag);
  Vario* variog           = Vario::computeFromDb(*varioparamG, grid, ECalcVario::VARIOGRAM);
  variog->display();
  delete variog;

  // With a quick entry
  mestitle(0, "Experimental variogram on Grid (quick entry)");
  auto* variog2 = varioGridCalculate(grid, ECalcVario::VARIOGRAM, true, 5);
  variog2->display();
  delete variog2;

  // ==========================================
  // Calculating Variogram Map on Isolated Data
  // ==========================================

  mestitle(0, "Variogram Map on Isolated Data");
  Db* vmapP = db_vmap(db);
  vmapP->display();
  delete vmapP;

  // =================================
  // Calculating Variogram Map on Grid
  // =================================

  mestitle(0, "Variogram Map on Grid");
  Db* vmapG = db_vmap(grid);
  DbStringFormat dbfmt(FLAG_STATS, {"VMAP*"});
  vmapG->display(&dbfmt);
  delete vmapG;

  delete db;
  delete grid;

  // ======================================
  // Calculation of all forms of Variograms
  // ======================================

  mestitle(0, "Calculation of all forms of Variograms on Standard Data");
  Id nvar     = 2;
  Id nsample  = 10;
  auto mesh   = 1.;
  auto* dbStd = DbGrid::createFillRandom({nsample}, nvar, 0, 0, 0., 0.,
                                         VectorDouble(), VectorDouble(), VectorDouble(), {mesh});
  dbStd->displayStats().display();

  auto nblag = 4;
  auto dlag  = mesh;
  for (const auto& name: ECalcVario::getAllKeys())
  {
    if (name == "UNDEFINED" || name == "GENERAL1" || name == "GENERAL2" || name == "GENERAL3") continue;
    for (bool flag_ergodic: {false, true})
    {
      auto* varioStd = variogramCalculate(dbStd, ECalcVario::fromKey(name), flag_ergodic, nblag, dlag);
      varioStd->display();
      delete varioStd;
    }
  }
  delete dbStd;

  // ===================================
  // Manipulation of variogram assessors
  // ===================================

  mestitle(0, "Manipulation of variogram assessors");

  // Create a 2-D data set and calculate an omni-directional variogram
  auto* dat        = Db::createFillRandom(10, 2, 2);
  auto* varioParam = VarioParam::createOmniDirection(10, 0.1);
  auto* vario      = Vario::computeFromDb(*varioParam, dat);
  vario->display();

  bool compress = true;

  // Extract the values from a multivariate variogram
  VectorDouble vec00 = vario->getGgVec(0, 0, 0, false, false, compress);
  printVector(vec00, "Variogram for: ivar=0, jvar=0 (no compress)", true, true);

  VectorDouble vec10 = vario->getGgVec(0, 1, 0, false, false, compress);
  printVector(vec10, "Variogram for: ivar=1, jvar=0 (no compress)", true, true);

  VectorDouble vec11 = vario->getGgVec(0, 1, 1, false, false, compress);
  printVector(vec11, "Variogram for: ivar=1, jvar=1 (no compress)", true, true);

  (void)vario->setGgVec(0, 0, 1, vec00);
  (void)vario->setGgVec(0, 1, 1, vec00);

  compress = false;

  vec00 = vario->getGgVec(0, 0, 0, false, false, compress);
  printVector(vec00, "Variogram for: ivar=0, jvar=0", true, true);

  vec10 = vario->getGgVec(0, 1, 0, false, false, compress);
  printVector(vec10, "Variogram for: ivar=1, jvar=0", true, true);

  vec11 = vario->getGgVec(0, 1, 1, false, false, compress);
  printVector(vec11, "Variogram for: ivar=1, jvar=1", true, true);

  (void)vario->setGgVec(0, 0, 1, vec00);
  (void)vario->setGgVec(0, 1, 1, vec00);
  vario->display();

  return static_cast<int>(error);
}
