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

// This test is mean to check the Factorial Kriging Analysis

#include "Basic/AStringFormat.hpp"
#include "Basic/NamingConvention.hpp"
#include "Covariances/CovContext.hpp"
#include "Enum/ESpaceType.hpp"

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCustom.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Estimation/CalcKriging.hpp"

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);
  AStringFormat format;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Parameters
  double oldstyle = 1.;
  bool verbose    = true;
  int nech        = 3;
  int nvar        = 2; // Should not be modified (see Model)
  bool flagSK     = true;
  OptCustom::define("oldStyle", oldstyle);

  // Generate the data base
  Db* data = Db::createFillRandom(nech, ndim, nvar, 0);
  data->setLocVariable(ELoc::Z, 1, 0, TEST);
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_ARRAY);
  data->display(dbfmt);

  // Generate the target file
  Db* target = Db::createFillRandom(1, ndim, 0, 0);
  target->setCoordinate(0, 0, data->getCoordinate(0, 0));
  target->setCoordinate(0, 1, data->getCoordinate(0, 1));

  // Create the Model
  Model* model = Model::create(CovContext(nvar, ndim));
  MatrixSquareSymmetric* sillNug =
    MatrixSquareSymmetric::createFromVD({1., 0.2, 0.2, 3.});
  model->addCovFromParam(ECov::NUGGET, 0., 0., 0., VectorDouble(), *sillNug);
  MatrixSquareSymmetric* sillExp =
    MatrixSquareSymmetric::createFromVD({2., 0.1, 0.1, 1.});
  model->addCovFromParam(ECov::EXPONENTIAL, 0.7, 0., 0., VectorDouble(), *sillExp);
  model->setCovFiltered(0, true);
  if (flagSK)
  {
    VectorDouble means = VH::simulateGaussian(nvar);
    model->setMeans(means);
  }
  else
    model->setDriftIRF(0, 0);
  VectorDouble means(nvar, 0.);

  // Neighborhood
  ANeigh* neigh = NeighUnique::create();

  // Define the verbose option
  if (verbose) OptDbg::setReference(1);

  // Test on Collocated CoKriging in Unique Neighborhood
  kriging(data, target, model, neigh);
  dbfmt = DbStringFormat::create(FLAG_STATS, {"Kriging.*"});
  target->display(dbfmt);

  // Free pointers

  delete sillNug;
  delete sillExp;
  delete neigh;
  delete data;
  delete target;
  delete model;

  return (0);
}
