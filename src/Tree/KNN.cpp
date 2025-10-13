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
#include "Tree/KNN.hpp"
#include "Basic/AStringable.hpp"
#include "Tree/ball_algorithm.h"

namespace gstlrn
{
KNN::KNN()
  : _heap()
  , _n_samples(0)
  , _n_neighbors(0)
{
}

bool KNN::_query(const t_btree& tree,
                 const MatrixT<double>& x,
                 Id n_samples,
                 Id n_features,
                 Id n_neigh)
{
  if (n_features != tree.n_features)
  {
    messerr(
      "query data dimension (%d) must match training data dimension (%d).",
      n_features, tree.n_features);
    _n_samples = -1;
    return false;
  }
  if (tree.n_samples < n_neigh)
  {
    messerr("'n_neigh' (%d) must be less than or equal to the number of "
            "training points (%d).",
            n_neigh, tree.n_samples);
    _n_samples = -1;
    return false;
  }
  _heap.resize(n_samples, n_neigh);
  _heap.load(tree, x);
  _heap.sort();

  return true;
}

Id KNN::btree_query(const t_btree& tree,
                     const MatrixT<double>& x,
                     Id n_samples,
                     Id n_features,
                     Id n_neigh)
{
  auto res = _query(tree, x, n_samples, n_features, n_neigh);
  if (!res) return 1;

  _n_samples   = _heap.n_pts;
  _n_neighbors = _heap.n_nbrs;

  return 0;
}

Id KNN::btree_query_inPlace(const t_btree& tree,
                             const MatrixT<double>& x,
                             Id n_samples,
                             Id n_features,
                             Id n_neigh,
                             Id rank,
                             VectorInt& indices,
                             VectorDouble& distances)
{
  if (rank < 0 || rank >= n_samples) return 1;

  auto res = _query(tree, x, n_samples, n_features, n_neigh);
  if (res)
  {
    indices   = {_heap.indices.getRow(rank).begin(), _heap.indices.getRow(rank).end()};
    distances = {_heap.distances.getRow(rank).begin(), _heap.distances.getRow(rank).end()};
  }
  return 0;
}

constvectint KNN::getIndices(Id rank) const
{
  if (rank < 0 || rank >= _n_samples) return {};
  return _heap.indices.getRow(rank);
}

Id KNN::getIndex(Id rank, Id ineigh) const
{
  if (rank < 0 || rank >= _n_samples) return ITEST;
  if (ineigh < 0 || ineigh >= _n_neighbors) return ITEST;
  return _heap.indices(rank, ineigh);
}

constvect KNN::getDistances(Id rank) const
{
  if (rank < 0 || rank >= _n_samples) return {};
  return _heap.distances.getRow(rank);
}

double KNN::getDistance(Id rank, Id ineigh) const
{
  if (rank < 0 || rank >= _n_samples) return ITEST;
  if (ineigh < 0 || ineigh >= _n_neighbors) return ITEST;
  return _heap.distances(rank, ineigh);
}
} // namespace gstlrn
