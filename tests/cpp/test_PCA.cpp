/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/*                                                                            */
/* This file is meant to demonstrate the PCA feature                          */
/*                                                                            */
/******************************************************************************/
#include "Enum/ECov.hpp"

#include "Space/ASpaceObject.hpp"
#include "Model/Model.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Stats/PCA.hpp"
#include "Stats/PCAStringFormat.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

#include <stdlib.h>

/****************************************************************************/
/*!
** Main Program for testing the sparse matrix algebra
**
*****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  double eps = 1.e-6;
  int error = 1;

  int ndim = 2;
  int nvar = 1;
  int nbsimu = 5;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  defineDefaultSpace(ESpaceType::RN, ndim);
  CovContext ctxt(nvar,ndim);

  // Creating a Point Data base in the 1x1 square with 'nech' samples
  int nech = 1000;
  Db* db = Db::createFromBox(nech,{0.,0.},{1.,1.}, 432423);

  // Creating the Model(s) of the Underlying GRF(s)
  double range1 = 0.2;
  Model* models = Model::createFromParam(ECov::BESSEL_K, range1, 1., 1.);
  models->display();

  // Perform a non-conditional simulation on the Db and on the Grid
  error = simtub(nullptr,db,models,nullptr,nbsimu);
  db->display();

  // ============
  // Evaluate PCA
  // ============

  mestitle(0,"Testing PCA");
  PCA pca = PCA(db);
  pca.display();

  // Store the transformed variables

  pca.dbZ2F(db);
  pca.dbF2Z(db);

  // Comparing initial and back-transformed variables
  VectorString names1 = generateMultipleNames("Simu" , nbsimu, ".");
  VectorString names2 = generateMultipleNames("F2Z.Z2F.Simu", nbsimu, ".");
  for (int i = 0; i < nbsimu; i++)
    (void) db->areSame(names1[i], names2[i], eps);

  // ============
  // Evaluate MAF
  // ============

  mestitle(0,"Testing MAF");
  db->setLocator("Simu*", ELoc::Z);
  DirParam dirparam = DirParam();
  PCA maf = PCA(db, 0.1, 0.05, dirparam);
  maf.display();

  // Store the transformed variable

  maf.dbZ2F(db, false, NamingConvention("Z2MAF"));
  maf.dbF2Z(db, false, NamingConvention("MAF2Z"));

  // Comparing initial and back-transformed variables

  names1 = generateMultipleNames("Simu" , nbsimu, ".");
  names2 = generateMultipleNames("MAF2Z.Z2MAF.Simu", nbsimu, ".");
  for (int i = 0; i < nbsimu; i++)
    (void) db->areSame(names1[i], names2[i], eps);

  delete db;
  return (error);
}
