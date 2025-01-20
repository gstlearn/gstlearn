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

// This test is mean to check the Collocated Cokriging
// - in Moving or Unique Neighborhood
// An additional printout is available

#include "Basic/AStringFormat.hpp"
#include "Basic/NamingConvention.hpp"
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
  double oldstyle = 0.;
  bool debug      = false;
  int nech        = 3;
  int nvar        = 3;
  bool flagSK     = false;
  bool flagUnique = true;
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
  Model* model = Model::createFromParam(ECov::EXPONENTIAL, scale, 0., 0., VectorDouble(),
                                 *sills, VectorDouble(), nullptr, false);
  if (flagSK)
  {
    VectorDouble means = VH::simulateGaussian(nvar);
    model->setMeans(means);
  }
  else
    model->setDriftIRF(0, 0);
  VectorDouble means(nvar, 0.);

  // Neighborhood
  ANeigh* neigh;
  int nmaxi     = nech;
  double radius = 5.;
  if (flagUnique)
    neigh = NeighUnique::create();
  else
    neigh = NeighMoving::create(false, nmaxi, radius);

  // Define the verbose option
  if (debug) OptDbg::setReference(1);

  // Test on Collocated CoKriging in Unique Neighborhood
  VectorInt varColCok = {0, 2};
  kriging(data, target, model, neigh, EKrigOpt::POINT, true, true, false,
          varColCok);
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
