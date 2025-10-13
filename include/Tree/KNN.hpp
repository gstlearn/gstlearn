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
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Tree/ball_algorithm.h"

namespace gstlrn
{

class GSTLEARN_EXPORT KNN
{
public:
  KNN();

  void setNNeighbors(Id n_neighbors) { _n_neighbors = n_neighbors; }
  void setNSamples(Id n_samples) { _n_samples = n_samples; }

  Id btree_query(const t_btree& tree,
                 const MatrixT<double>& x,
                 Id n_samples,
                 Id n_features,
                 Id n_neigh);
  Id btree_query_inPlace(const t_btree& tree,
                         const MatrixT<double>& x,
                         Id n_samples,
                         Id n_features,
                         Id n_neigh,
                         Id rank,
                         VectorInt& indices,
                         VectorDouble& distances);
#ifndef SWIG
  constvectint getIndices(Id rank = 0) const;
#endif // SWIG
  Id getIndex(Id rank = 0, Id ineigh = 0) const;
  constvect getDistances(Id rank = 0) const;
  double getDistance(Id rank = 0, Id ineigh = 0) const;

private:
  bool _query(const t_btree& tree,
              const MatrixT<double>& x,
              Id n_samples,
              Id n_features,
              Id n_neigh);

private:
  t_nheap _heap;
  Id _n_samples;
  Id _n_neighbors;
};
} // namespace gstlrn
