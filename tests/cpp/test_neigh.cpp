/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighWork.hpp"
#include "Basic/Vector.hpp"
#include "Basic/File.hpp"

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
  bool verbose = true;

  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  // Generate the data base
  int nech = 20;
  int seed = 43243;
  VectorDouble coormin = {0.,0.};
  VectorDouble coormax = {100.,100.};
  Db* db = Db::createFromBox(nech, coormin, coormax, ndim, seed, true);
  VectorDouble tab = ut_vector_simulate_gaussian(nech);
  db->addColumns(tab, "Variable", ELoc::Z);
  db->display();

  // Creating the target data base
  nech = 4;
  seed = 5436;
  Db* target = Db::createFromBox(nech, coormin, coormax, ndim, seed, true);
  target->display();

  // Creating a Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create(ndim, false);
  neighU->display();

  // Initializing the Neighborhood search
  mestitle(1,"Testing Unique Neighborhood");
  NeighWork nbghw(db,neighU);

  // Getting the Neighborhood for various target point
  ut_ivector_display("For Target Point #0",
                     nbghw.select(target, 0, VectorInt(), verbose));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());
  ut_ivector_display("For Target Point #1",
                     nbghw.select(target, 1, VectorInt(), verbose));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());
  delete neighU;

  // Creating a Moving Neighborhood
  int nmaxi = 5;
  double radius = 30.;
  NeighMoving* neighM = NeighMoving::create(ndim, false, nmaxi, radius);
  neighM->display();

  // Initializing the Neighborhood search
  mestitle(1,"Testing Moving Neighborhood");
  nbghw = NeighWork(db,neighM);

  // Getting the Neighborhood for various target point
  ut_ivector_display("For Target Point #0",
                     nbghw.select(target, 0, VectorInt(), verbose));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());
  ut_ivector_display("For Target Point #1",
                     nbghw.select(target, 1, VectorInt(), verbose));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());
  ut_ivector_display("For Target Point #2",
                     nbghw.select(target, 2, VectorInt(), verbose));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());
  ut_ivector_display("For Target Point #3",
                     nbghw.select(target, 3, VectorInt(), verbose));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());

  delete neighM;
  delete db;
  delete target;
  return (0);
}
