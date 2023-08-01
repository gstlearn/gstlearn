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

class Db;

class GSTLEARN_EXPORT Ball
{
public:
  Ball(double **data = nullptr,
       int n_samples = 0,
       int n_features = 0,
       int leaf_size = 10,
       int dist_type = 1);
  Ball(VectorVectorDouble &data, int leaf_size = 10, int dist_type = 1);
  Ball(const Db* db, int leaf_size = 10, int dist_type = 1, bool useSel = false);
  Ball(const Ball& p) = delete;
  Ball & operator=(const Ball& p) = delete;
  virtual ~Ball();

  t_knn query(double **test, int n_samples, int n_features, int n_neighbors = 1);
  t_knn queryAsVVD(VectorVectorDouble& test, int n_neighbors = 1);
  t_knn queryOne(double *test, int n_features, int n_neighbors = 1);
  t_knn queryOneAsVD(VectorDouble& test, int n_neighbors = 1);

protected:
  int _getFeatureNumber() const { return _tree->n_features; }
  int _getLeafSize() const { return _tree->leaf_size; }
  int _getSampleNumber() const { return _tree->n_samples; }
  int _getDistType() const { return _tree->dist_type; }
  bool _isValidFeatureNumber(int n_features) const;

private:
  mutable t_btree* _tree;
  mutable double (* _distFunction)(double, double , int);
};

GSTLEARN_EXPORT void  display(t_knn& knn, int ns_max = -1, int nn_max = -1);
GSTLEARN_EXPORT VectorInt getIndices(t_knn& knn, int rank = 0);
GSTLEARN_EXPORT VectorDouble getDistance(t_knn& knn, int rank = 0);
