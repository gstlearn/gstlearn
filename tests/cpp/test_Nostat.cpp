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
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"

#include <math.h>

#define __USE_MATH_DEFINES
#include <cmath>

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This test is meant to compare different manners to perform
 ** a non-conditional simulation using SPDE approach.
 ** The non-stationarity is provided through the angle of a Spiral which is povided:
 ** - directly as a pre-programmed non-stationary function (flagDirect = true)
 ** - by defining a vector containing the non-stationary angle (flagDirect = false and flagByAngle = true)
 ** - by defining the vector of local tensor matrices (flagDirect = false and flagByAngle = false)
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  int seed = 10355;
  law_set_random_seed(seed);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("test_nostat-");

  // Creating the 2-D Db
  auto nx = { 101, 101 };
  DbGrid* workingDbc = DbGrid::create(nx);

  // Creating the Non-stationary Model
  Model* model = Model::createFromParam(ECov::MATERN, 1., 1., 1., {10., 45.});

  FunctionalSpirale spirale(0., -1.4, 1., 1., 50., 50.);
  CovAniso* cova = model->getCovAniso(0);

  bool flagDirect  = true;
  bool flagByAngle = false;
  bool flagInquiry = false;

  if (flagDirect)
  {
    cova->makeAngleNoStatFunctional(&spirale);
  }
  else
  {
    if (flagByAngle)
    {
      VectorDouble angle = spirale.getFunctionValues(workingDbc);
      workingDbc->addColumns(angle, "Angle");
      cova->makeAngleNoStatDb("Angle",0,workingDbc);
    }
    else
    {
      VectorVectorDouble hh = spirale.getFunctionVectors(workingDbc, cova);
      workingDbc->addColumns(hh[0], "H1-1", ELoc::NOSTAT, 0);
      workingDbc->addColumns(hh[1], "H1-2", ELoc::NOSTAT, 1);
      workingDbc->addColumns(hh[2], "H2-2", ELoc::NOSTAT, 2);
      cova->makeTensorNoStatDb("H1-1",0,0,workingDbc);
      cova->makeTensorNoStatDb("H1-2",0,1,workingDbc);
      cova->makeTensorNoStatDb("H2-2",1,1,workingDbc);
    }
  }

  // Inquiry the value of the Non-stationary parameters at a given sample
  if (flagInquiry)
  {
    int target = 1000;
    VectorDouble vect = workingDbc->getSampleLocators(ELoc::NOSTAT, target);
    VH::dump("Non-stationary parameters at sample", vect);
  }

  int useCholesky = 0;
  law_set_random_seed(13256);
  (void)simulateSPDE(nullptr, workingDbc, model, 1, useCholesky,
                     VectorMeshes(), nullptr, VectorMeshes(), nullptr, SPDEParam(),
                     NamingConvention("Simu", true, false));

  DbStringFormat dbfmt(FLAG_STATS,{"Simu"});
  workingDbc->display(&dbfmt);

  (void) workingDbc->dumpToNF("spirale.ascii");

  message("Test performed successfully\n");

  delete workingDbc;
  delete model;
  return 0;
}
