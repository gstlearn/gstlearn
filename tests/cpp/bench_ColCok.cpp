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
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Enum/ESpaceType.hpp"
#include "Estimation/Estimations.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Space/ASpaceObject.hpp"

using namespace gstlrn;

void st_test(
  Id rank,
  bool flagSK,
  bool flagUnique,
  bool flagColCok,
  bool flagEnviron = false,
  Id iechref = 0,
  Id ndim = 2,
  Id nech = 3,
  Id nvar = 3)
{
  // Global parameters
  law_set_random_seed(32131);
  AStringFormat format;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Display parameters
  if (!flagEnviron)
  {
    message("\n Case %d", rank);
    message(" - %s", (flagSK) ? "Simple Kriging" : "Ordinary Kriging");
    message((flagUnique) ? " - Unique Neighborhood" : " - Moving Neighborhood");
    if (flagColCok) message(" - Collocation Option");
    message("\n");
  }

  // Define the Model
  Id order = (flagSK) ? -1 : 0;
  auto types = {ECov::EXPONENTIAL};
  Model* model = Model::createFillRandom(ndim, nvar, types, 1., order);
  VectorDouble means = VH::simulateGaussian(nvar);
  if (flagSK) model->setMeans(means);
  if (flagEnviron) model->display();

  // Define the Neighborhood
  ANeigh* neigh;
  if (flagUnique)
    neigh = NeighUnique::create();
  else
    neigh = NeighMoving::create(false, nech, 5.);
  if (flagEnviron) neigh->display();

  // Define the data base (with 'nech' samples)
  Db* data = Db::createFillRandom(nech, ndim, nvar, 0);
  data->setLocVariable(ELoc::Z, 1, 0, TEST);
  if (flagEnviron)
  {
    message("Data Db\n");
    data->getContentsAsTable().display();
  }

  // Define the target file (variables must also exist in Target for ColCok)
  // It contains two samples: the first is randomly located; the second coincides with the 'nech-1' data
  Db* target = Db::createFillRandom(
    2, ndim, nvar, 0, 0, 0., 0., VectorDouble(), VectorDouble(), VectorDouble(),
    534243);
  for (Id idim = 0; idim < ndim; idim++)
    target->setCoordinate(1, idim, data->getCoordinate(nech - 1, idim));
  if (flagEnviron)
  {
    message("Target\n");
    target->getContentsAsTable().display();
  }

  // Define the verbose option
  OptDbg::setReference(iechref);

  // Test on Collocated CoKriging in Unique Neighborhood
  VectorInt rank_colcok = {1, -1, 0};
  if (flagEnviron)
    rank_colcok.dump("Rank of ColCok variables (optional)", false);

  KrigOpt krigopt;
  if (flagColCok) krigopt.setColCok(rank_colcok);
  kriging(data, target, model, neigh, true, true, true, krigopt);

  // Print the results (only in the non verbose case, e.g. the active case)
  if (!flagEnviron) target->getContentsAsTable({"Kriging*"}).display();

  // Free pointers
  delete neigh;
  delete data;
  delete target;
  delete model;
}

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  st_test(0, true, false, false, true);

  // Perform the trials
  // st_test(rank, flagSK, flagUnique, flagColCok, flagEnviron);

  st_test(1, false, false, false);
  st_test(2, false, true, false);

  st_test(3, true, false, false);
  st_test(4, true, true, false);

  st_test(5, false, false, true);
  st_test(6, false, true, true);

  st_test(7, true, false, true);
  st_test(8, true, true, true);

  return (0);
}
