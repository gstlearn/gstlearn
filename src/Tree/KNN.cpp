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
#include "Tree/ball_algorithm.h"
#include "Basic/AStringable.hpp"

KNN::KNN()
  : _distances(),
    _indices(),
    _n_samples(0),
    _n_neighbors(0)
{
}

KNN::KNN(const KNN& m)
  : _distances()
  , _indices()
  , _n_samples(m._n_samples)
  , _n_neighbors(m._n_neighbors)
{
  _distances = m._distances;
  _indices   = m._indices;
}

KNN& KNN::operator=(const KNN &m)
{
  if (this != &m)
  {
    _n_samples   = m._n_samples;
    _n_neighbors = m._n_neighbors;
    _distances   = m._distances;
    _indices     = m._indices;
  }
  return *this;
}

KNN::~KNN()
{
}

t_nheap* KNN::_query(
  t_btree* tree, const double** x, int n_samples, int n_features, int n_neigh)
{
  t_nheap* heap = nullptr;

  if (n_features != tree->n_features)
  {
    messerr(
      "query data dimension (%d) must match training data dimension (%d).",
      n_features, tree->n_features);
    _n_samples = -1;
    return heap;
  }
  if (tree->n_samples < n_neigh)
  {
    messerr("'n_neigh' (%d) must be less than or equal to the number of "
            "training points (%d).",
            n_neigh, tree->n_samples);
    _n_samples = -1;
    return heap;
  }
  heap = nheap_init(n_samples, n_neigh);
  nheap_load(heap, tree, x);
  nheap_sort(heap);

  return heap;
}

int KNN::btree_query(
  t_btree* tree, const double** x, int n_samples, int n_features, int n_neigh)
{
  t_nheap* heap = _query(tree, x, n_samples, n_features, n_neigh);
  if (heap == nullptr) return 1;

  // Returned arguments
  _distances = copy_double_toVVD((const double**)heap->distances, heap->n_pts,
                                 heap->n_nbrs);
  _indices =
    copy_int_toVVI((const int**)heap->indices, heap->n_pts, heap->n_nbrs);
  _n_samples   = heap->n_pts;
  _n_neighbors = heap->n_nbrs;
  free(heap);

  return 0;
}

int KNN::btree_query_inPlace(t_btree* tree,
                             const double** x,
                             int n_samples,
                             int n_features,
                             int n_neigh,
                             int rank,
                             VectorInt& indices,
                             VectorDouble& distances)
{
  if (rank < 0 || rank >= n_samples) return 1;

  t_nheap* heap = _query(tree, x, n_samples, n_features, n_neigh);
  if (heap != nullptr)
  {
    int number = heap->n_nbrs;
    indices.resize(number);
    distances.resize(number);
    for (int j = 0; j < number; j++)
    {
      indices[j]   = heap->indices[rank][j];
      distances[j] = heap->distances[rank][j];
    }
    free(heap);
  }
  return 0;
}

VectorInt KNN::getIndices(int rank) const
{
  if (rank < 0 || rank >= _n_samples) return VectorInt();
  return _indices[rank];
}

int KNN::getIndex(int rank, int ineigh) const
{
  if (rank < 0 || rank >= _n_samples) return ITEST;
  if (ineigh < 0 || ineigh >= _n_neighbors) return ITEST;
  return _indices[rank][ineigh];
}

VectorDouble KNN::getDistances(int rank) const
{
  if (rank < 0 || rank >= _n_samples) return VectorDouble();
  return _distances[rank];
}

double KNN::getDistance(int rank, int ineigh) const
{
  if (rank < 0 || rank >= _n_samples) return ITEST;
  if (ineigh < 0 || ineigh >= _n_neighbors) return ITEST;
  return _distances[rank][ineigh];
}