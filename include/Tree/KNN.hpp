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

#include "Tree/ball_algorithm.h"
#include "Basic/VectorNumT.hpp"

struct t_btree;

class GSTLEARN_EXPORT KNN
{
public:
  KNN();
  KNN(const KNN& m);
  KNN & operator=(const KNN& m);
  virtual ~KNN();

  void setDistances(double **distances) { _distances = distances; }
  void setIndices(int **indices) { _indices = indices; }
  void setNNeighbors(int n_neighbors) { _n_neighbors = n_neighbors; }
  void setNSamples(int n_samples) { _n_samples = n_samples; }

  int btree_query(t_btree* tree,
                  const double **x,
                  int n_samples,
                  int n_features,
                  int n_neigh);
  VectorInt getIndices(int rank = 0) const;
  int getIndex(int rank = 0, int ineigh = 0) const;
  VectorDouble getDistance(int rank = 0) const;

private:
  double** _distances;
  int**    _indices;
  int      _n_samples;
  int      _n_neighbors;
};
