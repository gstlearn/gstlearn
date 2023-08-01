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
#include "Basic/AStringable.hpp"

Ball::Ball(double **data, int n_samples, int n_features, int leaf_size, int dist_type)
    : _tree(nullptr),
      _distFunction()
{
  _tree = btree_init(data, n_samples, n_features, leaf_size, dist_type);
}

Ball::Ball(VectorVectorDouble& data, int leaf_size, int dist_type)
    : _tree(nullptr),
      _distFunction()
{
  int n_samples = (int) data[0].size();
  int n_features = (int) data.size();
  _tree = btree_init(copy_double_arrAsVVD(data), n_samples, n_features, leaf_size, dist_type);
}

Ball::Ball(const Db* db, int leaf_size, int dist_type, bool useSel)
: _tree(nullptr),
  _distFunction()
{
  VectorVectorDouble data = db->getAllCoordinates(useSel);
  int n_samples = (int) data[0].size();
  int n_features = (int) data.size();
  _tree = btree_init(copy_double_arrAsVVD(data), n_samples, n_features, leaf_size, dist_type);
}

Ball::~Ball()
{
  free_tree(_tree);
}

t_knn Ball::query(double **test, int n_samples, int n_features, int n_neighbors)
{
  t_knn knn;
  if (! _isValidFeatureNumber(n_features)) return knn;
  knn = btree_query(_tree, test, n_samples, n_features, n_neighbors);
  return knn;
}

t_knn Ball::queryAsVVD(VectorVectorDouble& test, int n_neighbors)
{
  t_knn knn;
  if (test.empty()) return knn;
  int n_samples = (int) test[0].size();
  int n_features = (int) test.size();
  if (! _isValidFeatureNumber(n_features)) return knn;
  knn = btree_query(_tree, copy_double_arrAsVVD(test), n_samples, n_features, n_neighbors);
  return knn;
}

t_knn Ball::queryOne(double *test, int n_features, int n_neighbors)
{
  t_knn knn;
  if (! _isValidFeatureNumber(n_features)) return knn;
  double** internal = &test;
  knn = btree_query(_tree, internal, 1, n_features, n_neighbors);
  return knn;
}

t_knn Ball::queryOneAsVD(VectorDouble& test, int n_neighbors)
{
  t_knn knn;
  int n_features = (int) test.size();
  if (! _isValidFeatureNumber(n_features)) return knn;
  double* local = test.data();
  double **internal = &local;
  knn = btree_query(_tree, internal, 1, n_features, n_neighbors);
  return knn;
}

bool Ball::_isValidFeatureNumber(int n_features) const
{
  if (n_features != _getFeatureNumber())
  {
    messerr("Tree (%d) and Test (%d) must share the same feature dimension",
            n_features, _getFeatureNumber());
    return false;
  }
  return true;
}

void display(t_knn& knn, int ns_max, int nn_max)
{
  int ns = knn.n_samples;
  if (ns_max >= 0) ns = MIN(ns, ns_max);
  int nn = knn.n_neighbors;
  if (nn_max > 0) nn = MIN(nn, nn_max);
  for (int is = 0; is < ns; is++)
  {
    message("Indices = ");
    for (int in = 0; in < nn; in++)
      message(" %d", knn.indices[is][in]);
    message("\n");

    message("Distances = ");
    for (int in = 0; in < nn; in++)
      message(" %lf", knn.distances[is][in]);
    message("\n");
  }
}

VectorInt getIndices(t_knn& knn, int rank)
{
  if (rank < 0 || rank >= knn.n_samples) return VectorInt();
  return VectorHelper::initVInt(knn.indices[rank], knn.n_neighbors);
}
VectorDouble getDistance(t_knn& knn, int rank)
{
  if (rank < 0 || rank >= knn.n_samples) return VectorDouble();
  return VectorHelper::initVDouble(knn.distances[rank], knn.n_neighbors);
}
