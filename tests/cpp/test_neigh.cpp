/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/File.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program demonstrates the capabilities of Neigh classes
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Global parameters
  int ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  // Generate the data base
  int nech = 20;
  VectorDouble coormin = {0.,0.};
  VectorDouble coormax = {100.,100.};
  Db* db = Db::createFromBox(nech, coormin, coormax, 12345);
  VectorDouble tab = VH::simulateGaussian(nech);
  db->addColumns(tab, "Variable", ELoc::Z);
  db->display();

  // Creating the target data base
  nech = 4;
  Db* target = Db::createFromBox(nech, coormin, coormax, 12345);
  target->display();

  // Creating a Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create();
  neighU->attach(db, target);
  neighU->display();

  // Initializing the Neighborhood search
  mestitle(1,"Testing Unique Neighborhood");

  // Getting the Neighborhood for various target point
  VH::display("For Target Point #0", neighU->select(0));
  message("Is neighborhood Unchanged since last call = %d\n", neighU->isUnchanged());
  VH::display("For Target Point #1", neighU->select(1));
  message("Is neighborhood Unchanged since last call = %d\n", neighU->isUnchanged());
  delete neighU;

  // Creating a Moving Neighborhood
  int nmaxi = 5;
  double radius = 30.;
  NeighMoving* neighM = NeighMoving::create(false, nmaxi, radius);
  neighM->attach(db, target);
  neighM->display();

  // Initializing the Neighborhood search
  mestitle(1,"Testing Moving Neighborhood");

  // Getting the Neighborhood for various target point
  VH::display("For Target Point #0", neighM->select(0));
  message("Is neighborhood Unchanged since last call = %d\n", neighM->isUnchanged());
  VH::display("For Target Point #1", neighM->select(1));
  message("Is neighborhood Unchanged since last call = %d\n", neighM->isUnchanged());
  VH::display("For Target Point #2", neighM->select(2));
  message("Is neighborhood Unchanged since last call = %d\n", neighM->isUnchanged());
  VH::display("For Target Point #3", neighM->select(3));
  message("Is neighborhood Unchanged since last call = %d\n", neighM->isUnchanged());
  delete neighM;

  delete db;
  delete target;
  return (0);
}
