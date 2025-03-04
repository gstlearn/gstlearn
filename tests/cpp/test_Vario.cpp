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
#include "Covariances/ACov.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Enum/ECalcVario.hpp"
#include "Enum/ECov.hpp"

#include "Variogram/Vario.hpp"
#include "Variogram/VMap.hpp"
#include "Model/Model.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Covariances/CovAniso.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  int error = 1;
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  CovContext ctxt(1,2,1.); // use default space

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 1000;
  Db* db = Db::createFromBox(nech,{0.,0.},{1.,1.}, 3242);
  db->display();

  // Creating a grid covering the same space
  VectorInt nx = { 100, 100 };
  VectorDouble dx = { 0.01, 0.01 };
  DbGrid* grid = DbGrid::create(nx, dx);
  grid->display();

  // Creating the Model(s) of the Underlying GRF(s)
  Model models(ctxt);
  CovAnisoList covs(ctxt);
  double range1 = 0.2;
  CovAniso cova1(ECov::MATERN,range1,1.,1.,ctxt);
  covs.addCov(&cova1);
  models.setCovAnisoList(&covs);
  models.display();

  // Perform a non-conditional simulation on the Db and on the Grid
  error = simtub(nullptr,db,&models);
  db->display();
  error = simtub(nullptr,grid,&models);
  grid->display();

  // ===============
  // On Data samples
  // ===============

  mestitle(1, "Experimental variogram on Data Samples");
  int nlag = 20;
  VarioParam* varioparamP = VarioParam::createMultiple(2, nlag, 0.5 / nlag);
  Vario* variop = Vario::computeFromDb(*varioparamP,db,ECalcVario::VARIOGRAM);
  variop->display();
  message("Maximum Variogram Value = %lf\n",variop->getGmax());

  // Fitting the experimental variogram of Underlying GRF (with constraint that total sill is 1)
  Model model(ctxt);
  VectorECov covas {ECov::MATERN, ECov::EXPONENTIAL};
  model.fit(variop,covas,false);
  model.display();

  // ===============
  // On Grid samples
  // ===============

  mestitle(1, "Experimental variogram on Grid");
  VarioParam* varioparamG = VarioParam::createMultipleFromGrid(grid, nlag);
  Vario* variog = Vario::computeFromDb(*varioparamG, grid, ECalcVario::VARIOGRAM);
  variog->display();

  // ==========================================
  // Calculating Variogram Map on Isolated Data
  // ==========================================

  mestitle(1, "Variogram Map on Isolated Data");
  Db* vmapP = db_vmap(db, ECalcVario::VARIOGRAM);
  vmapP->display();

  // =================================
  // Calculating Variogram Map on Grid
  // =================================

  mestitle(1, "Variogram Map on Grid");
  Db* vmapG = db_vmap(grid, ECalcVario::VARIOGRAM);
  DbStringFormat dbfmt(FLAG_STATS,{"VMAP*"});
  vmapG->display(&dbfmt);

  delete db;
  delete grid;
  delete varioparamP;
  delete variop;
  delete variog;
  delete vmapG;
  delete vmapP;

  return (error);
}
