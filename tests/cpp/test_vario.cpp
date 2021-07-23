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
#include "Neigh/Neigh.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"
#include "Basic/AStringable.hpp"

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int argc, char *argv[])

{
  auto pygst = std::string(std::getenv("PYGSTLEARN_DIR"));
  int error = 0;
  int ndim = 2;
  CovContext ctxt(1,2,1.);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 1000;
  Db db(nech,VectorDouble(2,0.),VectorDouble(2,1.));
  db.display(1);

  // Creating a grid covering the same space
  VectorInt nx = { 100, 100 };
  VectorDouble dx = { 0.01, 0.01 };
  Db grid(nx, dx);
  grid.display(1);

  // Creating the Model(s) of the Underlying GRF(s)
  Model models(ctxt);
  double range1 = 0.2;
  CovAniso cova1(COV_BESSEL_K,range1,1.,1.,ctxt);
  models.addCova(&cova1);
  models.display();

  // Perform a non-conditional simulation on the Db and on the Grid

  error = simtub(nullptr,&db,&models);
  db.display(1);
  error = simtub(nullptr,&grid,&models);
  grid.display();

  // ===============
  // On Data samples
  // ===============

  // Determination of the experimental variogram
  Vario variod;
  int nlag = 20;
  std::vector<Dir> dirs = generateMultipleDirs(ndim, 2, nlag, 0.5 / nlag);
  variod.addDirs(dirs);
  error = variod.compute(&db,"vg");
  variod.display(1);
  message("Maximum Variogram Value = %lf\n",variod.getGmax());

  // Fitting the experimental variogram o Underlying GRF (with constraint that total sill is 1)
  Model model(ctxt);
  std::vector<ENUM_COVS> covs {COV_BESSEL_K, COV_EXPONENTIAL};
  model.fit(&variod,covs,true);
  model.display();

  // ===============
  // On Grid samples
  // ===============

  // Determination of the experimental variogram
  Vario variog;
  std::vector<Dir> dirgs = generateMultipleGridDirs(ndim, nlag);
  variog.addDirs(dirgs);
  error = variog.compute(&grid,"vg",VectorDouble(),VectorDouble(),true);
  variog.display(1);

  return(error);
}
