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

// This test is meant to check the Cross-validation estimation
// - in mutlivariate case
// - in Moving and Unique Neighborhood
// A printout is provided

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/NamingConvention.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCst.hpp"
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
  bool flagSK     = false;
  bool flagUnique = false;
  bool verbose    = true;
  int nech        = 3;
  int nvar        = (flagUnique) ? 1 : 2;

  // Generate the data base
  Db* data = Db::createFillRandom(nech, ndim, nvar, 0);
  data->setLocVariable(ELoc::Z, 1, 0, TEST);
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_ARRAY);
  data->display(dbfmt);

  // Create the Model
  double scale = 0.7;
  MatrixSquareSymmetric* sills = MatrixSquareSymmetric::createRandomDefinitePositive(nvar);
  Model* model =
    Model::createFromParam(ECov::EXPONENTIAL, scale, 0., 0., VectorDouble(),
                           *sills, VectorDouble(), nullptr, false);
  if (flagSK)
  {
    VectorDouble means = VH::simulateGaussian(nvar);
    model->setMeans(means);
  }
  else
    model->setDriftIRF(0, 0);

  // Unique Neighborhood
  ANeigh* neigh;
  int nmaxi     = nech;
  double radius = 5.;
  if (flagUnique)
    neigh = NeighUnique::create();
  else
    neigh = NeighMoving::create(false, nmaxi, radius);

  // Define the verbose option
  if (verbose) OptDbg::setReference(1);
  OptCst::define(ECst::NTCOL, -1);

  // Perform the cross-validation
  data->setLocator("z*", ELoc::Z);
  xvalid(data, model, neigh, false);
  OptDbg::setReference(0);

  // Perform the cross-validation
  data->setLocator("z*", ELoc::Z);
  xvalid(data, model, neigh, false, -1, -1);

  // Produce some statistics for comparison
  dbfmt = DbStringFormat::create(FLAG_ARRAY, {"Xvalid.*"});
  data->display(dbfmt);

  // Free pointers

  delete sills;
  delete neigh;
  delete data;
  delete model;
  delete dbfmt;

  return (0);
}
