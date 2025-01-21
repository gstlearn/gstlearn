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

// This test is meant to check CoKriging with Bayesian hypotheses on Drift
// - in Multivariate (heterotopic) case
// - in Unique or Moving neighborhood

#include "Basic/NamingConvention.hpp"
#include "Enum/ESpaceType.hpp"

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Space/ASpace.hpp"
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
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Parameters
  double oldstyle = 0.;
  bool verbose    = false;
  int nech        = 3;
  int nvar        = 1;
  OptCustom::define("oldStyle", oldstyle);

  // Generate the data base
  Db* data = Db::createFillRandom(nech, ndim, nvar, 0);
  data->setLocVariable(ELoc::Z, 1, 0, TEST);
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_ARRAY);
  data->display(dbfmt);

  // Generate the target file
  Db* target = Db::createFillRandom(1, ndim, 0, 0);

  // Create the Model
  double scale = 0.7;
  MatrixSquareSymmetric* sills =
    MatrixSquareSymmetric::createRandomDefinitePositive(nvar);
  Model* model =
    Model::createFromParam(ECov::EXPONENTIAL, scale, 0., 0., VectorDouble(),
                           *sills, VectorDouble(), ASpaceSharedPtr(), false);
  model->setDriftIRF(0, 0);

  // Neighborhood
  ANeigh* neigh = NeighUnique::create();

  // Create the Bayesian Priors for Drift coefficients
  VectorDouble PriorMean = VH::simulateGaussian(nvar);
  MatrixSquareSymmetric PriorCov(nvar);
  PriorCov.setDiagonal(VH::simulateUniform(nvar, 0.1, 0.5));

  // Define the verbose option
  if (verbose) OptDbg::setReference(1);

  // Test on Bayesian
  kribayes(data, target, model, neigh, PriorMean, PriorCov);
  dbfmt = DbStringFormat::create(FLAG_STATS, {"Kriging.*"});
  target->display(dbfmt);

  // Free pointers

  delete sills;
  delete neigh;
  delete data;
  delete target;
  delete model;

  return (0);
}
