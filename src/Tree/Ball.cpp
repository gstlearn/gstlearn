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
#include "Tree/ball_algorithm.h"
#include "Db/Db.hpp"
#include "Space/SpacePoint.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"

Ball::Ball(const double** data,
           int n_samples,
           int n_features,
           double (*dist_function)(const double* x1,
                                   const double* x2,
                                   int size),
           int leaf_size,
           int default_distance_function)
  : _tree(nullptr)
{
  _tree = btree_init(data, n_samples, n_features, false, dist_function, leaf_size,
                     default_distance_function);
}

Ball::Ball(const VectorVectorDouble& data,
           double (*dist_function)(const double* x1,
                                   const double* x2,
                                   int size),
           int leaf_size,
           int default_distance_function)
  : _tree(nullptr)
{
  int n_samples     = (int)data[0].size();
  int n_features    = (int)data.size();
  double** internal = copy_double_arrAsVVD(data);
  _tree = btree_init((const double**)internal, n_samples, n_features, false,
                     dist_function, leaf_size, default_distance_function);
  free_2d_double(internal, n_features);
}

Ball::Ball(const Db* dbin,
           const Db* dbout,
           double (*dist_function)(const double* x1,
                                   const double* x2,
                                   int size),
           int leaf_size,
           bool has_constraints,
           int default_distance_function,
           bool useSel)
  : _tree(nullptr)
{
 
  VectorVectorDouble data = dbin->getAllCoordinates(useSel);
  double** internal       = copy_double_arrAsVVD(data);
  int n_samples           = (int)data[0].size();
  int n_features          = (int)data.size();

  if (dbout != nullptr)
  {
    VectorVectorDouble aux = dbout->getAllCoordinates(useSel);
    append_double_arrAsVVD(aux, &internal, n_samples);
    n_samples += (int)aux[0].size();
  }

  _tree = btree_init((const double**)internal, n_samples, n_features, has_constraints,
                     dist_function, leaf_size, default_distance_function);
  free_2d_double(internal, n_samples);
}

void Ball::init(const Db* db,
                double (*dist_function)(const double* x1,
                                        const double* x2,
                                        int size),
                int leaf_size,
                int default_distance_function,
                bool useSel)
{
  if (_tree != nullptr) free_tree(_tree);
  VectorVectorDouble data = db->getAllCoordinates(useSel);
  int n_samples           = (int)data[0].size();
  int n_features          = (int)data.size();
  double** internal       = copy_double_arrAsVVD(data);
  _tree = btree_init((const double**)internal, n_samples, n_features, false,
                     dist_function, leaf_size, default_distance_function);
  // free_2d_double(internal, n_features);
  free_2d_double(internal, n_samples);
}

Ball::~Ball()
{
  free_tree(_tree);
}

KNN Ball::query(const double **test, int n_samples, int n_features, int n_neighbors)
{
  KNN knn;
  (void) knn.btree_query(_tree, test, n_samples, n_features, n_neighbors);
  return knn;
}

KNN Ball::queryAsVVD(const VectorVectorDouble& test, int n_neighbors)
{
  KNN knn;
  if (test.empty()) return knn;
  int n_samples = (int) test[0].size();
  int n_features = (int) test.size();
  double** internal = copy_double_arrAsVVD(test);
  (void) knn.btree_query(_tree, (const double**) internal, n_samples, n_features, n_neighbors);
  //  free_2d_double(internal, n_features);
  free_2d_double(internal, n_samples);
  return knn;
}

KNN Ball::queryOne(const double *test, int n_features, int n_neighbors)
{
  KNN knn;
  (void) knn.btree_query(_tree, (const double**) &test, 1, n_features, n_neighbors);
  return knn;
}

KNN Ball::queryOneAsVD(const VectorDouble& test, int n_neighbors)
{
  KNN knn;
  int n_features = (int) test.size();
  const double* internal = test.data();
  (void) knn.btree_query(_tree, (const double**) &internal, 1, n_features, n_neighbors);
  return knn;
}

KNN Ball::queryOneAsVDFromSP(const SpacePoint& Pt, int n_neighbors)
{
  KNN knn;
  int n_features = Pt.getNDim();
  const double* internal = Pt.getCoords().data();
  (void)knn.btree_query(_tree, (const double**)&internal, 1, n_features,
                        n_neighbors);
  return knn;
}

VectorInt Ball::getIndices(const SpacePoint& Pt, int n_neighbors)
{
  KNN knn = queryOneAsVDFromSP(Pt, n_neighbors);
  return knn.getIndices(0);
}

int Ball::queryClosest(const VectorDouble& test)
{
  KNN knn;
  int n_features = (int) test.size();
  const double* internal = test.data();
  if (knn.btree_query(_tree, (const double**) &internal, 1, n_features, 1)) return ITEST;
  return knn.getIndex(0, 0);
}

int Ball::queryOneInPlace(const VectorDouble& test,
                          int n_neighbors,
                          VectorInt& indices,
                          VectorDouble& distances,
                          int rank)
{
  KNN knn;
  int n_features         = (int)test.size();
  const double* internal = test.data();
  return knn.btree_query_inPlace(_tree, (const double**)&internal, 1,
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
void Ball::display(int level) const
{
  btree_display(_tree, level);
}

bool Ball::_isConstraintDefined() const
{
  if (_tree->accept == nullptr)
  {
    messerr("You may not set one Constraint if not initialized in Ball constructor");
    return false;
  }
  return true;
}

int Ball::setConstraint(int rank, bool status)
{
  if (_tree == nullptr) return 1;
  if (! _isConstraintDefined()) return 1;
  if (rank < 0 || rank >= _tree->n_samples) return 1;
  _tree->accept[rank] = status;
  return 0;
}

int Ball::resetConstraints(bool status)
{
  if (_tree == nullptr) return 1;
  if (!_isConstraintDefined()) return 1;
  for (int i = 0, n = _tree->n_samples; i < n; i++)
    _tree->accept[i] = status;
  return 0;
}

MatrixT<int> findNN(Db* dbin,
                    Db* dbout,
                    int nb_neigh,
                    bool flagShuffle,
                    bool verbose,
                    double (*dist_function)(const double* x1,
                                            const double* x2,
                                            int size),
                    int leaf_size,
                    int default_distance_function)
{
  MatrixT<int> mat;

  // Preliminary checks
  int ndim = dbin->getNDim();
  if (dbout != nullptr && ndim != dbout->getNDim())
  {
    messerr("Dbin(%d) and Dbout(%d) should have the same dimension",
            ndim, dbout->getNDim());
    return mat;
  }

  // Creating the Ball tree
  Ball ball(dbin, dbout, dist_function, leaf_size, true, default_distance_function);
  if (verbose) ball.display(1);

  // Dimensioning the output matrix
  int n1 = dbin->getNSample(true);
  int n2 = (dbout != nullptr) ? dbout->getNSample(true) : 0;
  mat.resize(n1 + n2, nb_neigh);

  // Loop on the samples for the FNN search
  SpacePoint pt;
  VectorInt ranks;
  VectorInt neighs(nb_neigh);
  VectorDouble distances(nb_neigh);

  if (verbose)
    mestitle(1, "List of Neighborhoors for NN search");

  ranks = (flagShuffle) ? law_random_path(n1) : VH::sequence(n1);
  for (int jech = 0; jech < n1; jech++)
  {
    int iech = ranks[jech];
    dbin->getSampleAsSPInPlace(pt, iech);
    ball.setConstraint(iech, true);
    (void)ball.queryOneInPlace(pt.getCoordUnprotected(), nb_neigh, neighs, distances);
    for (int i = 0; i < nb_neigh; i++) mat(jech, i) = neighs[i];

    if (verbose)
    {
      message("Sample_1 %3d", iech);
      VH::dump(" ", neighs, false);
    }
  }

  if (dbout != nullptr)
  {
    ranks = (flagShuffle) ? law_random_path(n2) : VH::sequence(n2);
    for (int jech = 0; jech < n2; jech++)
    {
      int iech = ranks[jech];
      dbout->getSampleAsSPInPlace(pt, iech);
      ball.setConstraint(iech + n1, true);
      (void)ball.queryOneInPlace(pt.getCoordUnprotected(), nb_neigh, neighs, distances);
      for (int i = 0; i < nb_neigh; i++) mat(n1 + jech, i) = neighs[i];

      if (verbose)
      {
        message("Sample_2 %3d", n1 + iech);
        VH::dump(" ", neighs, false);
      }
    }
  }
  return mat;
}
