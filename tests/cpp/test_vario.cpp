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
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Covariances/ECov.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Variogram/ECalcVario.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

#include <stdlib.h>

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  int error = 1;
  int ndim = 2;
  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, ndim);
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
  std::vector<DirParam> dirparamP = DirParam::createMultiple(ndim, 2, nlag, 0.5 / nlag);
  varioparamP.addMultiDirs(dirparamP);
  Vario variop = Vario(&varioparamP,db);
  variop.computeByKey("vg");
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
  std::vector<DirParam> dirparamG = DirParam::createMultipleFromGrid(ndim, nlag);
  varioparamG.addMultiDirs(dirparamG);
  Vario variog = Vario(&varioparamG, grid);
  variog.computeByKey("vg");
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
  DbStringFormat dbfmt(FLAG_STATS,{"VMAP*"});
  vmapG->display(&dbfmt);

  delete db;
  delete grid;
  return (error);
}
