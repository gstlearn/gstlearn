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

Ball::Ball(const double **data, int n_samples, int n_features, int leaf_size, int dist_type)
    : _tree(nullptr)
{
  _tree = btree_init(data, n_samples, n_features, leaf_size, dist_type);
}

Ball::Ball(const VectorVectorDouble& data, int leaf_size, int dist_type)
    : _tree(nullptr)
{
  int n_samples = (int) data[0].size();
  int n_features = (int) data.size();
  double** internal = copy_double_arrAsVVD(data);
  _tree = btree_init((const double**) internal, n_samples, n_features,
                     leaf_size, dist_type);
  free_2d_double(internal, n_features);
}

Ball::Ball(const Db *db, int leaf_size, int dist_type, bool useSel)
    : _tree(nullptr)
{
  VectorVectorDouble data = db->getAllCoordinates(useSel);
  int n_samples = (int) data[0].size();
  int n_features = (int) data.size();
  double** internal = copy_double_arrAsVVD(data);
  _tree = btree_init((const double**) internal, n_samples,
                     n_features, leaf_size, dist_type);
  free_2d_double(internal, n_features);
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
  free_2d_double(internal, n_features);
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

int Ball::queryClosest(const VectorDouble& test)
{
  KNN knn;
  int n_features = (int) test.size();
  const double* internal = test.data();
  if (knn.btree_query(_tree, (const double**) &internal, 1, n_features, 1)) return ITEST;
  return knn.getIndex(0, 0);
}

void Ball::display(int level) const
{
  btree_display(_tree, level);
}
