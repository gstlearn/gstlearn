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

VectorDouble
getSortedDistance(Db* data, const VectorInt& ranks, const SpaceTarget& Pt)
{
  int size = (int)ranks.size();
  SpaceTarget Dt;
  VectorDouble dist(size, TEST);
  for (int iech = 0; iech < size; iech++)
  {
    data->getSampleAsSTInPlace(ranks[iech], Dt);
    dist = Pt.getDistance(Dt);
  }
  VH::sortInPlace(dist);
  return dist;
}

/****************************************************************************/
/*!
 ** Main Program
 ** This is meant to compare the Moving Neighborhood to the KNN features
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("NeighKNN-");

  // Global parameters
  defineDefaultSpace(ESpaceType::RN, 2);
  Timer timer;

  // Main parameters
  int ndat      = 100000;
  int ntarget   = 1000;
  int nmaxi     = 20;
  int leaf_size = 20;
  int mode      = 0;
  bool verbose = false;

  // Generate the input data base
  int ndim = 2;
  int nvar = 1;
  Db* data = Db::createFillRandom(ndat, ndim, nvar);

  // Generate the output data base
  Db* target  = Db::createFillRandom(ntarget, ndim, 0);

  // Moving Neighborhood
  double radius = 0.5;
  NeighMoving* neigh1 = NeighMoving::create(false, nmaxi, radius);
  neigh1->attach(data, target);

  NeighMoving* neigh2 = NeighMoving::create(false, nmaxi, radius);
  neigh2->setBallSearch(true, leaf_size);
  neigh2->attach(data, target);

  // General printout
  mestitle(1, "Neighborhood search");
  message("Comparing Moving Neighborhood searches\n");
  message("- with Standard search\n");
  message("- with BallTree search\n");
  message("General characteristics\n");
  message("- Number of data = %d\n", ndat);
  message("- Number of targets = %d\n", ntarget);
  message("- Number of neighboring samples = %d\n", nmaxi);
  message("- Leaf Size for KNN = %d\n", leaf_size);
  message("Computing time includes Data to Target distance evaluation for comparison\n");

  SpaceTarget Pt;
  VectorVectorDouble checkDistances1;
  VectorVectorDouble checkDistances2;

  // Using Standard neighborhood
  if (mode <= 0 || mode == 1)
  {
    message("\nStandard Neighborhood search\n");

    timer.reset();
    if (verbose) message("Data Indices for each Target Neighborhood\n");
    VectorInt indices1;
    for (int i = 0; i < ntarget; i++)
    {
      target->getSampleAsSTInPlace(i, Pt);
      neigh1->getNeigh(i, indices1);
      VectorDouble dists = getSortedDistance(data, indices1, Pt);
      checkDistances1.push_back(dists);
      if (verbose) VH::dump("", indices1);
    }
    timer.displayIntervalMilliseconds("Standard Neighborhood", 17000);
  }

  // Using Neighborhood with Ball Tree option
  if (mode <= 0 || mode == 2)
  {
    message("\nNeighborhood search with Ball-Tree option\n");

    timer.reset();
    if (verbose) message("Data Indices for each Target Neighborhood\n");
    VectorInt indices2;
    for (int i = 0; i < ntarget; i++)
    {
      target->getSampleAsSTInPlace(i, Pt);
      neigh2->getNeigh(i, indices2);
      VectorDouble dists = getSortedDistance(data, indices2, Pt);
      checkDistances2.push_back(dists);
      if (verbose) VH::dump("", indices2);
    }
    timer.displayIntervalMilliseconds("Neigh. with Ball option", 400);
  }

  // Compare the set of indices
  if (mode <= 0)
  {
    for (int i = 0; i < ntarget; i++)
      if (!VH::isEqual(checkDistances1[i], checkDistances2[i], 0.01))
      {
        messerr("Vector of indices are different at rank %d", i);
        VH::dump("- Standard search", checkDistances1[i]);
        VH::dump("- BallTree search", checkDistances2[i]);
      }
  }

  delete neigh1;
  delete neigh2;
  delete data;
  delete target;

  return (0);
  }
