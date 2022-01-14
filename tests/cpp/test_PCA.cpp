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
/* This file is meant to demonstrate the PCA feature                 */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"
#include "Model/Model.hpp"
#include "Basic/AStringable.hpp"
#include "Db/Db.hpp"
#include "Covariances/ECov.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Stats/PCA.hpp"
#include <stdlib.h>

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  int error = 1;

  int ndim = 2;
  int nvar = 1;
  int nbsimu = 5;

  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  CovContext ctxt(nvar,ndim);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 1000;
  Db* db = Db::createFromBox(nech,VectorDouble(2,0.),VectorDouble(2,1.));

  // Creating the Model(s) of the Underlying GRF(s)
  Model models(ctxt);
  CovLMC covs(ctxt.getSpace());
  double range1 = 0.2;
  CovAniso cova1(ECov::BESSEL_K,range1,1.,1.,ctxt);
  covs.addCov(&cova1);
  models.setCovList(&covs);
  models.display();

  // Perform a non-conditional simulation on the Db and on the Grid
  error = simtub(nullptr,db,&models,nullptr,nbsimu);
  db->display();

  // ===============
  // Evaluate PCA
  // ===============

  PCA pca(nbsimu);
  pca.compute(db);
  pca.display();

  // Store the transformed variables

  pca.dbZ2F(db);
  db->display();

  delete db;
  return (error);
}
