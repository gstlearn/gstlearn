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
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/File.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighWork.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program demonstrates the capabilities of Neigh and NeighWork classes
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
  neighU->display();

  // Initializing the Neighborhood search
  mestitle(1,"Testing Unique Neighborhood");
  NeighWork* neighW = NeighWork::create(db,neighU,target);

  // Getting the Neighborhood for various target point
  VH::display("For Target Point #0", neighW->select(0));
  message("Is neighborhood Unchanged since last call = %d\n", neighW->isUnchanged());
  VH::display("For Target Point #1", neighW->select(1));
  message("Is neighborhood Unchanged since last call = %d\n", neighW->isUnchanged());
  delete neighU;
  delete neighW;

  // Creating a Moving Neighborhood
  int nmaxi = 5;
  double radius = 30.;
  NeighMoving* neighM = NeighMoving::create(false, nmaxi, radius);
  neighM->display();

  // Initializing the Neighborhood search
  mestitle(1,"Testing Moving Neighborhood");
  neighW = NeighWork::create(db,neighM,target);

  // Getting the Neighborhood for various target point
  VH::display("For Target Point #0", neighW->select(0));
  message("Is neighborhood Unchanged since last call = %d\n", neighW->isUnchanged());
  VH::display("For Target Point #1", neighW->select(1));
  message("Is neighborhood Unchanged since last call = %d\n", neighW->isUnchanged());
  VH::display("For Target Point #2", neighW->select(2));
  message("Is neighborhood Unchanged since last call = %d\n", neighW->isUnchanged());
  VH::display("For Target Point #3", neighW->select(3));
  message("Is neighborhood Unchanged since last call = %d\n", neighW->isUnchanged());
  delete neighM;
  delete neighW;

  delete db;
  delete target;
  return (0);
}
