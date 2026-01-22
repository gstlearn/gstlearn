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
#include "Basic/File.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"

using namespace gstlrn;
/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program demonstrates the capabilities of Neigh classes
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  // Global parameters
  Id ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Generate the data base
  Id nech              = 20;
  VectorDouble coormin = {0., 0.};
  VectorDouble coormax = {100., 100.};
  Db* db               = Db::createFromBox(nech, coormin, coormax, 12345);
  VectorDouble tab     = VH::simulateGaussian(nech);
  db->addColumns(tab, "Variable", ELoc::Z);
  db->display();

  // Creating the target data base
  nech       = 4;
  Db* target = Db::createFromBox(nech, coormin, coormax, 12345);
  target->display();

  // Creating a Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create();
  neighU->attach(db, target);
  neighU->display();
  VectorInt nbgh;

  // Initializing the Neighborhood search
  mestitle(1, "Testing Unique Neighborhood");

  // Getting the Neighborhood for various target point
  nbgh.clear();
  neighU->select(0, nbgh);
  printVector(nbgh, "For Target Point #0", true, true);
  message("Is neighborhood Unchanged since last call = %d\n", neighU->isUnchanged());
  neighU->select(1, nbgh);
  printVector(nbgh, "For Target Point #1", true, true);
  message("Is neighborhood Unchanged since last call = %d\n", neighU->isUnchanged());
  delete neighU;

  // Creating a Moving Neighborhood
  Id nmaxi            = 5;
  double radius       = 30.;
  NeighMoving* neighM = NeighMoving::create(false, nmaxi, radius);
  neighM->attach(db, target);
  neighM->display();

  // Initializing the Neighborhood search
  mestitle(1, "Testing Moving Neighborhood");

  // Getting the Neighborhood for various target point
  nbgh.clear();
  neighM->select(0, nbgh);
  printVector(nbgh, "For Target Point #0", true, true);
  message("Is neighborhood Unchanged since last call = %d\n", neighM->isUnchanged());
  neighM->select(1, nbgh);
  printVector(nbgh, "For Target Point #1", true, true);
  message("Is neighborhood Unchanged since last call = %d\n", neighM->isUnchanged());
  neighM->select(2, nbgh);
  printVector(nbgh, "For Target Point #2", true, true);
  message("Is neighborhood Unchanged since last call = %d\n", neighM->isUnchanged());
  neighM->select(3, nbgh);
  printVector(nbgh, "For Target Point #3", true, true);
  message("Is neighborhood Unchanged since last call = %d\n", neighM->isUnchanged());
  delete neighM;

  delete db;
  delete target;
  return (0);
}
