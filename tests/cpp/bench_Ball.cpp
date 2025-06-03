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

#include "LinearOp/ProjMatrix.hpp"
#include "Mesh/MeshSphericalExt.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpaceSN.hpp"
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

  ASerializable::setPrefixName("Tree-");

  // Global parameters
  bool verbose = false;
  int ndim = 2;
  int mode = 0;
  defineDefaultSpace(ESpaceType::RN, ndim);

  // Constructing the Data Set
  int nech = 100;
  Db* data = Db::createFillRandom(nech, ndim, 1, 0, 0, 0., 0.,
                                  VectorDouble(), VectorDouble(),
                                  VectorDouble(), 131343);
  if (verbose) data->display();

  int nb_neigh = 5;
  VectorInt neighs;
  VectorDouble distances;

  if (mode == 0 || mode == 1)
  {
    // ================
    // Traditional Ball
    // ================
    mestitle(0, "Traditional use of the Ball Tree");

    // Constructing the Ball Tree
    Ball ball(data, nullptr, nullptr, 10, false);
    if (verbose) ball.display(0);

    // My target sample
    VectorDouble target = {0.4, 0.2};
    SpacePoint pt1(target);

    // Inquiring the Ball tree
    mestitle(1, "Various ways of inquiring the Ball Tree");

    // - for the closest sample
    int ineigh = ball.queryClosest(target);
    message("The closest sample to the Target is : %d\n", ineigh);

    // - for a set of neighboring samples
    neighs = ball.getIndices(pt1, nb_neigh);
    VH::dump("Indices of the target neighbors", neighs);

    // - for a more complete output (in place)
    (void)ball.queryOneInPlace(target, nb_neigh, neighs, distances);
    VH::dump("Indices of the neighbors", neighs);
    VH::dump("Distances to the target", distances);
  }

  if (mode == 0 || mode == 2)
  {
    // =====================
    // Ball with constraints
    // =====================
    mestitle(0, "Use of the Ball Tree with Constraints (FNN search)");
    bool has_constraints = true;
    verbose              = true;

    // Constructing the Ball Tree from Db(s)
    Ball ball(data, nullptr, nullptr, 10, has_constraints);
    if (verbose) ball.display(1);

    // Loop on the samples for the FNN search
    SpacePoint pt2;
    // VectorInt ranks = law_random_path(nech);
    VectorInt ranks = VH::sequence(nech);
    for (int jech = 0; jech < nech; jech++)
    {
      int iech = ranks[jech];
      data->getSampleAsSPInPlace(pt2, iech);
      ball.setConstraint(iech, true);
      (void)ball.queryOneInPlace(pt2.getCoordUnprotected(), nb_neigh, neighs, distances);
      VH::dump("Indices of the neighbors", neighs);
    }
  }

  if (mode == 0 || mode == 3)
  {
    // ================
    // FindNN algorithm
    // ================
    mestitle(0, "Demonstrating the findNN algorithm");
    bool flagShuffle = true;

    int nech = 20;
    Db* aux = Db::createFillRandom(nech, ndim, 1, 0, 0, 0., 0.,
                                    VectorDouble(), VectorDouble(),
                                    VectorDouble(), 24813);
    MatrixT<int> mat = findNN(data, aux, nb_neigh, flagShuffle);
    int nrows        = mat.getNRows();
    int ncols        = mat.getNCols();
    for (int irow = 0; irow < nrows; irow++)
    {
      for (int icol = 0; icol < ncols; icol++)
      {
        int value = mat(irow, icol);
        if (IFFFF(value))
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
    bool verbose = false;
    defineDefaultSpace(ESpaceType::SN, 2);

    // Constructing the Meshing on the Sphere
    MeshSphericalExt mesh = MeshSphericalExt();
    mesh.resetFromDb(NULL, NULL, "-r3");
    if (verbose) mesh.display();

    // Define the Ball Tree starting from the mesh
    Ball ball(&mesh);
    if (verbose) ball.display(1);

    // Defining a Target point in Longitude / Latitude
    VectorDouble target = {70., 70.};
    SpacePoint pt1(target);

    // Inquiring the Ball tree
    mestitle(1, "Various ways of inquiring the Ball Tree");

    // - for the closest sample
    int ineigh = ball.queryClosest(target);
    message("The closest sample to the Target is : %d\n", ineigh);

    // - for a set of neighboring samples
    (void)ball.queryOneInPlace(target, nb_neigh, neighs, distances);
    VH::dump("Indices of the neighbors", neighs);
    VH::dump("Distances to the target", distances);
  }
  return (0);
}
