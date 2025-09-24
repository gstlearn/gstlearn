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
#include "Tree/Ball.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Mesh/AMesh.hpp"
#include "Space/SpacePoint.hpp"
#include "Tree/ball_algorithm.h"

namespace gstlrn
{
/**
 * @brief Construct a new Ball tree object
 *
 * @param dbin  Input database
 * @param dbout Output database (appended to input) (can be null)
 * @param leaf_size Number of elements in the leafs of the Ball tree
 * @param all_available True if samples are available for selection at the beginning (True by default)
 * @param default_distance_function 1 for Euclidean distance, 2 for Manhattan
 * @param useSel True if only the selected samples of dbin are used to build the Ball tree
 */
Ball::Ball(const Db* dbin,
           const Db* dbout,
           Id leaf_size,
           bool all_available,
           Id default_distance_function,
           bool useSel)
  : _tree()
{
  Id n_samples;
  Id n_features;
  auto internal = _getInformationFromDb(dbin, dbout, useSel, &n_samples, &n_features);
  if (internal.empty()) return;

  _tree = t_btree(std::move(internal), n_samples, n_features, all_available,
                  leaf_size, default_distance_function);
}

/**
 * @brief Construct a new Ball tree object based on the barycenters of the meshes
 *
 * @param mesh  AMesh description
 * @param leaf_size Number of elements in the leafs of the Ball tree
 * @param all_available True if samples are available for selection at the beginning
 * @param default_distance_function 1 for Euclidean distance, 2 for Manhattan
 */
Ball::Ball(const AMesh* mesh,
           Id leaf_size,
           bool all_available,
           Id default_distance_function)
{
  Id n_samples;
  Id n_features;
  auto internal = _getInformationFromMesh(mesh, &n_samples, &n_features);
  if (internal.empty()) return;

  _tree = t_btree(std::move(internal), n_samples, n_features, all_available,
                  leaf_size, default_distance_function);
}
/**
 * @brief Construct a new Ball tree object
 *
 * @param db  Input database
 * @param leaf_size Number of elements in the leafs of the Ball tree
 * @param all_available True if samples are available for selection at the beginning (True by default)
 * @param default_distance_function 1 for Euclidean distance, 2 for Manhattan
 * @param useSel True if only the selected samples of dbin are used to build the Ball tree
 */
void Ball::init(const Db* db,
                Id leaf_size,
                bool all_available,
                Id default_distance_function,
                bool useSel)
{
  Id n_samples;
  Id n_features;
  auto internal = _getInformationFromDb(db, nullptr, useSel, &n_samples, &n_features);
  if (internal.empty()) return;

  _tree = t_btree(std::move(internal), n_samples, n_features, all_available, leaf_size, default_distance_function);
}

KNN Ball::queryAsVVD(const VectorVectorDouble& test, Id n_neighbors)
{
  KNN knn;
  if (test.empty()) return knn;
  Id n_samples  = static_cast<Id>(test[0].size());
  Id n_features = static_cast<Id>(test.size());
  MatrixT<double> internal(n_samples, n_features);
  // transpose
  for (Id i = 0; i < n_samples; ++i)
  {
    for (Id j = 0; j < n_features; ++j)
    {
      internal(i, j) = test[j][i];
    }
  }
  (void)knn.btree_query(_tree, internal, n_samples, n_features, n_neighbors);
  return knn;
}

KNN Ball::queryOne(const double* test, Id n_features, Id n_neighbors)
{
  KNN knn;
  MatrixT<double> internal(1, n_features);
  for (Id i = 0; i < n_features; ++i)
  {
    internal(0, i) = test[i];
  }
  (void)knn.btree_query(_tree, internal, 1, n_features, n_neighbors);
  return knn;
}

KNN Ball::queryOneAsVD(const VectorDouble& test, Id n_neighbors)
{
  KNN knn;
  Id n_features = static_cast<Id>(test.size());
  MatrixT<double> internal(1, n_features);
  for (Id i = 0; i < n_features; ++i)
  {
    internal(0, i) = test[i];
  }
  (void)knn.btree_query(_tree, internal, 1, n_features, n_neighbors);
  return knn;
}

KNN Ball::queryOneAsVDFromSP(const SpacePoint& Pt, Id n_neighbors)
{
  KNN knn;
  Id n_features = static_cast<Id>(Pt.getNDim());
  MatrixT<double> internal(1, n_features);
  for (Id i = 0; i < n_features; ++i)
  {
    internal(0, i) = Pt.getCoord(i);
  }
  (void)knn.btree_query(_tree, internal, 1, n_features, n_neighbors);
  return knn;
}

VectorInt Ball::getIndices(const SpacePoint& Pt, Id n_neighbors)
{
  KNN knn = queryOneAsVDFromSP(Pt, n_neighbors);
  return {knn.getIndices(0).begin(), knn.getIndices(0).end()};
}

Id Ball::queryClosest(const VectorDouble& test)
{
  KNN knn;
  Id n_features = static_cast<Id>(test.size());
  MatrixT<double> internal(1, n_features);
  for (Id i = 0; i < n_features; ++i)
  {
    internal(0, i) = test[i];
  }
  if (knn.btree_query(_tree, internal, 1, n_features, 1)) return ITEST;
  return knn.getIndex(0, 0);
}

Id Ball::queryOneInPlace(const VectorDouble& test,
                         Id n_neighbors,
                         VectorInt& indices,
                         VectorDouble& distances,
                         Id rank)
{
  KNN knn;
  Id n_features = static_cast<Id>(test.size());
  MatrixT<double> internal(1, n_features);
  for (Id i = 0; i < n_features; ++i)
  {
    internal(0, i) = test[i];
  }
  return knn.btree_query_inPlace(_tree, internal, 1,
                                 n_features, n_neighbors, rank, indices, distances);
}

/**
 * @brief Ask for information regarding the Ball Tree organization
 *
 * @param level Level of details
 *              -1 Just the general volumetry information
 *               0 List of the different nodes
 *               1 List of Leaves and attached list of samples
 */
void Ball::display(Id level) const
{
  _tree.display(level);
}

Id Ball::setAvailable(Id rank, bool status)
{
  if (empty()) return 1;
  if (rank < 0 || rank >= _tree.n_samples) return 1;
  _tree.available[rank] = status;
  return 0;
}

Id Ball::resetAvailable(bool status)
{
  if (empty()) return 1;
  for (Id i = 0, n = _tree.n_samples; i < n; i++)
    _tree.available[i] = status;
  return 0;
}

MatrixT<Id> findNN(const Db* dbin,
                   const Db* dbout,
                   Id nb_neigh,
                   bool flagShuffle,
                   bool verbose,
                   Id leaf_size,
                   Id default_distance_function)
{
  MatrixT<Id> mat;

  // Preliminary checks
  Id ndim = dbin->getNDim();
  if (dbout != nullptr && ndim != dbout->getNDim())
  {
    messerr("Dbin(%d) and Dbout(%d) should have the same dimension",
            ndim, dbout->getNDim());
    return mat;
  }

  // Creating the Ball tree
  Ball ball(dbin, dbout, leaf_size, false, default_distance_function);
  if (verbose) ball.display(1);

  // Dimensioning the output matrix
  Id n1 = dbin->getNSample(true);
  Id n2 = (dbout != nullptr) ? dbout->getNSample(true) : 0;
  mat.resize(n1 + n2, nb_neigh);

  // Loop on the samples for the FNN search
  SpacePoint pt;
  VectorInt ranks;
  VectorInt neighs(nb_neigh);
  VectorDouble distances(nb_neigh);

  if (verbose)
    mestitle(1, "List of Neighbors for NN search");

  ranks = (flagShuffle) ? law_random_path(n1) : VH::sequence(n1);
  for (Id jech = 0; jech < n1; jech++)
  {
    Id iech = ranks[jech];
    dbin->getSampleAsSPInPlace(pt, iech);
    ball.setAvailable(iech, true);
    (void)ball.queryOneInPlace(pt.getCoordsUnprotected(), nb_neigh, neighs, distances);
    for (Id i = 0; i < nb_neigh; i++) mat(jech, i) = neighs[i];

    if (verbose)
    {
      message("Sample_1 %3d", iech);
      VH::dump(" ", neighs, false);
    }
  }

  if (dbout != nullptr)
  {
    ranks = (flagShuffle) ? law_random_path(n2) : VH::sequence(n2);
    for (Id jech = 0; jech < n2; jech++)
    {
      Id iech = ranks[jech];
      dbout->getSampleAsSPInPlace(pt, iech);
      ball.setAvailable(iech + n1, true);
      (void)ball.queryOneInPlace(pt.getCoordsUnprotected(), nb_neigh, neighs, distances);
      for (Id i = 0; i < nb_neigh; i++) mat(n1 + jech, i) = neighs[i];

      if (verbose)
      {
        message("Sample_2 %3d", n1 + iech);
        VH::dump(" ", neighs, false);
      }
    }
  }
  return mat;
}

MatrixT<double> Ball::_getInformationFromDb(const Db* dbin,
                                            const Db* dbout,
                                            bool useSel,
                                            Id* n_samples,
                                            Id* n_features)
{
  VectorDouble oneColumn;
  Id ncol    = dbin->getNLoc(ELoc::X);
  Id nrowtot = dbin->getNSample(useSel);
  if (dbout != nullptr)
  {
    if (ncol != dbout->getNLoc(ELoc::X))
    {
      messerr("'dbin' and 'dbout' should share the same space dimension");
      return {};
    }
    nrowtot += dbout->getNSample(useSel);
  }

  // Core allocation
  MatrixT<double> internal(nrowtot, ncol);

  // Loading the information from dbin
  Id ns = 0;
  if (dbin != nullptr)
  {
    Id nrow = dbin->getNSample(useSel);
    for (Id icol = 0; icol < ncol; icol++)
    {
      oneColumn = dbin->getOneCoordinate(icol, useSel);
      for (Id irow = 0; irow < nrow; irow++)
        internal(ns + irow, icol) = oneColumn[irow];
    }
    ns += nrow;
  }

  // Loading information from dbout
  if (dbout != nullptr)
  {
    Id nrow = dbout->getNSample(useSel);
    for (Id icol = 0; icol < ncol; icol++)
    {
      oneColumn = dbout->getOneCoordinate(icol, useSel);
      for (Id irow = 0; irow < nrow; irow++)
        internal(ns + irow, icol) = oneColumn[irow];
    }
    ns += nrow;
  }

  *n_samples  = ns;
  *n_features = ncol;
  return internal;
}

MatrixT<double> Ball::_getInformationFromMesh(const AMesh* mesh,
                                              Id* n_samples,
                                              Id* n_features)
{
  VectorDouble oneColumn;
  Id ndim  = mesh->getNDim();
  Id nmesh = mesh->getNMeshes();

  // Core allocation
  MatrixT<double> internal(nmesh, ndim);

  // Loading the information from mesh
  for (Id imesh = 0; imesh < nmesh; imesh++)
  {
    mesh->getBarycenterInPlace(imesh, internal.getRow(imesh));
  }

  *n_samples  = nmesh;
  *n_features = ndim;
  return internal;
}
} // namespace gstlrn
