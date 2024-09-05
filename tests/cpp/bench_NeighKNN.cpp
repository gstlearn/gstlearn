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
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Basic/VectorHelper.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Tree/Ball.hpp"
#include "Tree/KNN.hpp"

/****************************************************************************/
/*!
 ** Main Program
 ** This is meant to compare the Moving Neighborhood to the KNN features
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  bool verbose = false;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("NeighKNN-");

  // Global parameters
  defineDefaultSpace(ESpaceType::RN, 2);
  Timer timer;

  // Main parameters
  int ndat = 100000;
  int ntarget = 1000;
  int nmaxi = 20;

  // Generate the input data base
  int ndim = 2;
  int nvar = 1;
  Db* data = Db::createFillRandom(ndat, ndim, nvar);

  // Generate the output data base
  Db* target  = Db::createFillRandom(ntarget, ndim, 0);

  // Moving Neighborhood
  double radius = 0.5;
  NeighMoving* neigh = NeighMoving::create(false, nmaxi, radius);
  neigh->attach(data, target);

  // General printout
  mestitle(1, "Neighborhood search");
  message("- Number of data = %d\n", ndat);
  message("- Number of targets = %d\n", ntarget);
  message("- Number of neighboring samples = %d\n", nmaxi);
  
  // Get the neighborhood
  message("\nStandard Neighborhood search\n");

  timer.reset();
  if (verbose) message("Data Indices for each Target Neighborhood\n");
  VectorInt ranks;
  for (int i = 0; i < ntarget; i++)
  {
    neigh->getNeigh(i, ranks);
    if (verbose)
    {
      VH::sortInPlace(ranks);
      VH::display("", ranks);
    }
  }
  timer.displayIntervalMilliseconds("Standard Neighborhood", 19000);

  // Using KNN
  message("\nKNN Neighborhood search\n");

  timer.reset();
  if (verbose) message("Data Indices for each Target Neighborhood\n");
  int leaf_size = 30;
  Ball ball(data, leaf_size);
  VectorDouble coor(ndim);
  VectorInt indices;
  for (int i = 0; i < ntarget; i++)
  {
    target->getCoordinatesPerSampleInPlace(i, coor);
    KNN knn = ball.queryOneAsVD(coor, nmaxi);
    indices = knn.getIndices(0);
    if (verbose)
    {
      VH::sortInPlace(indices);
      VH::display("", indices);
    }
  }
  timer.displayIntervalMilliseconds("KNN Neighborhood", 900);

  delete neigh;
  delete data;
  delete target;

  return (0);
  }
