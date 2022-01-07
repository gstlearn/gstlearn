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
#include "Neigh/Neigh.hpp"
#include "Neigh/NeighWork.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program demsontrates the capabilities of Neigh and NeighWork classes
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])

{
  // Global parameters
  int ndim = 2;
  int nvar = 1;

  // Generate the data base
  int nech = 20;
  int seed = 43243;
  VectorDouble coormin = {0.,0.};
  VectorDouble coormax = {100.,100.};
  Db db(nech, coormin, coormax, ndim, seed, true);
  VectorDouble tab = ut_vector_simulate_gaussian(nech);
  db.addFields(tab, "Variable", ELoc::Z);
  db.display();

  // Creating the target data base
  nech = 4;
  seed = 5436;
  Db target(nech, coormin, coormax, ndim, seed, true);
  target.display();

  // Creating a Unique Neighborhood
  Neigh neigh(ndim);
  neigh.display();

  // Initializing the Neighborhood search
  mestitle(1,"Testing Unique Neighborhood");
  NeighWork nbghw(&db,&neigh);

  // Getting the Neighborhood for various target point
  ut_ivector_display("For Target Point #0",
                     nbghw.select(&target, 0));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());
  ut_ivector_display("For Target Point #1",
                     nbghw.select(&target, 1));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());

  // Creating a Moving Neighborhood
  int nmaxi = 5;
  double radius = 30.;
  neigh = Neigh(ndim, nmaxi, radius);
  neigh.display();

  // Initializing the Neighborhood search
  mestitle(1,"Testing Moving Neighborhood");
  nbghw = NeighWork(&db,&neigh);

  // Getting the Neighborhood for various target point
  ut_ivector_display("For Target Point #0",
                     nbghw.select(&target, 0));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());
  ut_ivector_display("For Target Point #1",
                     nbghw.select(&target, 1));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());
  ut_ivector_display("For Target Point #2",
                     nbghw.select(&target, 2));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());
  ut_ivector_display("For Target Point #3",
                     nbghw.select(&target, 3));
  message("Is neighborhood Unchanged since last call = %d\n",
          nbghw.isUnchanged());


  return (0);
}
