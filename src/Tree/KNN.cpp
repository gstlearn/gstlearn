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
#include "Basic/VectorHelper.hpp"
#include "Basic/AStringable.hpp"

KNN::KNN()
  : _distances(nullptr),
    _indices(nullptr),
    _n_samples(0),
    _n_neighbors(0)
{
}

KNN::KNN(const KNN &m)
  : _distances(nullptr),
    _indices(nullptr),
    _n_samples(m._n_samples),
    _n_neighbors(m._n_neighbors)
{
  _distances = copy_double_arr((const double**) m._distances, m._n_samples, m._n_neighbors);
  _indices = copy_int_arr((const int**) m._indices, m._n_samples, m._n_neighbors);
}

KNN& KNN::operator=(const KNN &m)
{
  if (this != &m)
  {
    _n_samples = m._n_samples;
    _n_neighbors = m._n_neighbors;
    _distances = copy_double_arr((const double**) m._distances, m._n_samples, m._n_neighbors);
    _indices = copy_int_arr((const int**) m._indices, m._n_samples, m._n_neighbors);
  }
  return *this;
}

KNN::~KNN()
{
  free_2d_double(_distances, _n_samples);
  free_2d_int(_indices, _n_samples);
}

int KNN::btree_query(t_btree* tree, const double **x, int n_samples, int n_features, int n_neigh)
{
  t_nheap *heap;

  if (n_features != tree->n_features)
  {
    messerr("query data dimension (%d) must match training data dimension (%d).",
            n_features, tree->n_features);
    _n_samples = -1;
    return 1;
  }
  if (tree->n_samples < n_neigh)
  {
    messerr("'n_neigh' (%d) must be less than or equal to the number of training points (%d).",
            n_neigh, tree->n_samples);
    _n_samples = -1;
    return 1;
  }
  heap = nheap_init(n_samples, n_neigh);
  nheap_load(heap, tree, x);
  nheap_sort(heap);

  // Returned arguments
  _distances = copy_double_arr((const double **) heap->distances, heap->n_pts, heap->n_nbrs);
  _indices = copy_int_arr((const int **) heap->indices, heap->n_pts, heap->n_nbrs);
  _n_samples = heap->n_pts;
  _n_neighbors = heap->n_nbrs;
  free(heap);

  return 0;
}

VectorInt KNN::getIndices(int rank) const
{
  if (rank < 0 || rank >= _n_samples) return VectorInt();
  return VectorHelper::initVInt(_indices[rank], _n_neighbors);
}

int KNN::getIndex(int rank, int ineigh) const
{
  if (rank < 0 || rank >= _n_samples) return ITEST;
  if (ineigh < 0 || ineigh >= _n_neighbors) return ITEST;
  return _indices[rank][ineigh];
}

VectorDouble KNN::getDistance(int rank) const
{
  if (rank < 0 || rank >= _n_samples) return VectorDouble();
  return VectorHelper::initVDouble(_distances[rank], _n_neighbors);
}
