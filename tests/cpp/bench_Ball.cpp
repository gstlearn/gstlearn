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

#include "Basic/File.hpp"
#include "Basic/Timer.hpp"
#include "Basic/VectorHelper.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Db/Db.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "Mesh/MeshSphericalExt.hpp"
#include "Model/Model.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpaceSN.hpp"
#include "Tree/Ball.hpp"

using namespace gstlrn;
/****************************************************************************/
/*!
 ** Main Program
 ** This is meant to exhibit the Ball tree mechanism (for future improvements)
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  Timer timer;
  VectorDouble vec;
  VectorDouble vecb;
  VectorDouble diff;

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("bench_Ball-");

  // Global parameters
  bool verbose = false;
  Id ndim      = 2;
  Id mode      = 0;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Constructing the Data Set
  Id nech  = 100;
  Db* data = Db::createFillRandom(nech, ndim, 1, 0, 0, 0., 0.,
                                  VectorDouble(), VectorDouble(),
                                  VectorDouble(), 131343);
  if (verbose) data->display();

  Id nb_neigh = 5;
  VectorInt neighs;
  VectorDouble distances;

  if (mode == 0 || mode == 1)
  {
    // ================
    // Traditional Ball
    // ================
    mestitle(0, "Traditional use of the Ball Tree");

    // Constructing the Ball Tree
    Ball ball(data, nullptr, 10);
    if (verbose) ball.display(0);

    // My target sample
    VectorDouble target = {0.4, 0.2};
    SpacePoint pt1(target); // Use default space

    // Inquiring the Ball tree
    mestitle(1, "Various ways of inquiring the Ball Tree");

    // - for the closest sample
    Id ineigh = ball.queryClosest(target);
    message("The closest sample to the Target is : %d\n", ineigh);

    // - for a set of neighboring samples
    neighs = ball.getIndices(pt1, nb_neigh);
    printVector(neighs, "Indices of the target neighbors", true, true);

    // - for a more complete output (in place)
    (void)ball.queryOneInPlace(target, nb_neigh, neighs, distances);
    printVector(neighs, "Indices of the neighbors", true, true);
    printVector(distances, "Distances to the target", true, true);
  }

  if (mode == 0 || mode == 2)
  {
    // ===========================
    // Ball with 'Available flags'
    // ===========================
    mestitle(0, "Use of the Ball Tree with Constraints (FNN search)");
    bool all_available = false;
    verbose            = true;

    // Constructing the Ball Tree from Db(s)
    Ball ball(data, nullptr, 10, all_available);
    if (verbose) ball.display(1);

    // Loop on the samples for the FNN search
    SpacePoint pt2; // Use default space
    // VectorInt ranks = law_random_path(nech);
    VectorInt ranks = VH::sequence(nech);
    for (Id jech = 0; jech < nech; jech++)
    {
      Id iech = ranks[jech];
      data->getSampleAsSPInPlace(pt2, iech);
      ball.setAvailable(iech, true);
      (void)ball.queryOneInPlace(pt2.getCoordsUnprotected(), nb_neigh, neighs, distances);
      printVector(neighs, "Indices of the neighbors", true, true);
    }
  }

  if (mode == 0 || mode == 3)
  {
    // ================
    // FindNN algorithm
    // ================
    mestitle(0, "Demonstrating the findNN algorithm");
    bool flagShuffle = true;

    nech    = 20;
    Db* aux = Db::createFillRandom(nech, ndim, 1, 0, 0, 0., 0.,
                                   VectorDouble(), VectorDouble(),
                                   VectorDouble(), 24813);

    auto mat   = findNN(data, aux, nb_neigh, flagShuffle);
    auto nrows = static_cast<Id>(mat.getNRows());
    auto ncols = static_cast<Id>(mat.getNCols());
    for (Id irow = 0; irow < nrows; irow++)
    {
      for (Id icol = 0; icol < ncols; icol++)
      {
        Id value = mat(irow, icol);
        if (isNA(value))
          message("   NA");
        else
          message(" %4d", mat(irow, icol));
      }
      message("\n");
    }
  }

  // Cleaning
  delete data;

  if (mode == 0 || mode == 4)
  {
    // ===================
    // Working on a Sphere
    // ===================
    mestitle(0, "Demonstrating BallTree algorithm on a Sphere");

    // Global parameters
    verbose = false;
    defineDefaultSpace(ESpaceType::SN, 2);

    // Constructing the Meshing on the Sphere
    MeshSphericalExt mesh;
    mesh.resetFromDb(NULL, NULL, "-r3");
    if (verbose) mesh.display();

    // Define the Ball Tree starting from the mesh
    Ball ball(&mesh);
    if (verbose) ball.display(1);

    // Defining a Target point in Longitude / Latitude
    VectorDouble target = {70., 70.};
    SpacePoint pt1(target); // Use default space

    // Inquiring the Ball tree
    mestitle(1, "Various ways of inquiring the Ball Tree");

    // - for the closest sample
    Id ineigh = ball.queryClosest(target);
    message("The closest sample to the Target is : %d\n", ineigh);

    // - for a set of neighboring samples
    (void)ball.queryOneInPlace(target, nb_neigh, neighs, distances);
    printVector(neighs, "Indices of the neighbors", true, true);
    printVector(distances, "Distances to the target", true, true);
  }
  return (0);
}
