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
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  int ndim = 2;
  int nvar = 1;
  law_set_random_seed(32131);

  defineDefaultSpace(ESpaceType::RN, ndim);

  // Generate the data base
  int nech = 10;
  Db* data = Db::createFillRandom(nech, ndim, nvar);
  data->display();

  // Generate the target file
  Db* target = Db::createFillRandom(1, ndim, 0);
  target->display();

  // Create the Model
  double range = 0.5;
  double sill = 2.;
  Model* model = Model::createFromParam(ECov::SPHERICAL, range, sill);
  model->display();

  // Unique Neighborhood
  NeighUnique* neigh = NeighUnique::create();
  neigh->display();

  // Define the verbose option
  OptDbg::setReference(1);

  // ====================== Testing Neighborhood Storage ===========================
  message("\n---> Using Standard Kriging procedure\n");
  kriging(data, target, model, neigh);

  // ====================== Free pointers ==================================
  delete neigh;
  delete data;
  delete target;
  delete model;

  return (0);
}
