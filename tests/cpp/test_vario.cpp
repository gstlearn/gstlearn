/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/*                                                                            */
/* This file is meant to demonstrate the process of using PGS                 */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"
#include "Variogram/Vario.hpp"
#include "Model/Model.hpp"
#include "Basic/AStringable.hpp"
#include "Db/Db.hpp"
#include "Covariances/ECov.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include <stdlib.h>

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  int error = 1; //TODO : temporary fail
  int ndim = 2;
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  CovContext ctxt(1,2,1.); // use default space

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 1000;
  Db* db = Db::createFromBox(nech,VectorDouble(2,0.),VectorDouble(2,1.));
  db->display();

  // Creating a grid covering the same space
  VectorInt nx = { 100, 100 };
  VectorDouble dx = { 0.01, 0.01 };
  Db* grid = Db::createFromGrid(nx, dx);
  grid->display();

  // Creating the Model(s) of the Underlying GRF(s)
  Model models(ctxt);
  CovLMC covs(ctxt.getSpace());
  double range1 = 0.2;
  CovAniso cova1(ECov::BESSEL_K,range1,1.,1.,ctxt);
  covs.addCov(&cova1);
  models.setCovList(&covs);
  models.display();

  // Perform a non-conditional simulation on the Db and on the Grid
  error = simtub(nullptr,db,&models);
  db->display();
  error = simtub(nullptr,grid,&models);
  grid->display();

  // ===============
  // On Data samples
  // ===============

  // Determination of the experimental variogram
  VarioParam varioparamP;
  int nlag = 20;
  std::vector<DirParam> dirparamP = generateMultipleDirs(ndim, 2, nlag, 0.5 / nlag);
  varioparamP.addMultiDirs(dirparamP);
  Vario variop = Vario(&varioparamP,db);
  variop.compute("vg");
  variop.display();
  message("Maximum Variogram Value = %lf\n",variop.getGmax());

  // Fitting the experimental variogram of Underlying GRF (with constraint that total sill is 1)
  Model model(ctxt);
  std::vector<ECov> covas {ECov::BESSEL_K, ECov::EXPONENTIAL};
  model.fit(&variop,covas,true);
  model.display();

  // ===============
  // On Grid samples
  // ===============

  // Determination of the experimental variogram
  VarioParam varioparamG;
  std::vector<DirParam> dirparamG = generateMultipleGridDirs(ndim, nlag);
  varioparamG.addMultiDirs(dirparamG);
  Vario variog = Vario(&varioparamG, grid);
  variog.compute("vg",true);
  variog.display();

  // ==========================================
  // Calculating Variogram Map on Isolated Data
  // ==========================================

  Db* vmapP = db_vmap_compute(db, ECalcVario::VARIOGRAM);
  vmapP->display();

  // =================================
  // Calculating Variogram Map on Grid
  // =================================

  Db* vmapG = db_vmap_compute(grid, ECalcVario::VARIOGRAM);
  vmapG->display();

  delete db;
  delete grid;
  return (error);
}
