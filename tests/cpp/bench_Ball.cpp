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
#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Basic/VectorHelper.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Tree/Ball.hpp"

/****************************************************************************/
/*!
 ** Main Program
 ** This is meant to exhibit the Ball tree mechanism (for future improvements)
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  Timer timer;
  VectorDouble vec;
  VectorDouble vecb;
  VectorDouble diff;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Tree-");

  // Global parameters
  int ndim = 2;
  bool verbose = true;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Constructing the Data Set
  int nech = 100;
  Db* data = Db::createFillRandom(nech, ndim, 1, 0, 0, 0., 0.,
                                  VectorDouble(), VectorDouble(),
                                  VectorDouble(), 131343);
  if (verbose) data->display();

  // Constructing the Ball Tree
  Ball ball(data);
  if (verbose) ball.display(0);

  // Inquiring the Ball tree
  mestitle(0, "Various ways of inquiring the Ball Tree");

  // - for the closest sample
  VectorDouble coor = {0.4, 0.2};
  int ineigh = ball.queryClosest(coor);
  message("The closest sample to the Target is : %d\n", ineigh);

  // - for a set of neighboring samples
  int nb_neigh = 5;
  SpacePoint pt(coor);
  VectorInt neighs = ball.getIndices(pt, nb_neigh);
  VH::dump("Indices of the target neighbors", neighs);

  // -for a more complete output (in place)
  VectorDouble distances;
  (void)ball.queryOneInPlace(coor, nb_neigh, neighs, distances);
  VH::dump("Indices of the target neighbors", neighs);
  VH::dump("Distances to the target", distances);

  // Cleaning
  delete data;

  return (0);
}
