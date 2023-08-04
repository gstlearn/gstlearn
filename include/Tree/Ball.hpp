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
  Ball(const double **data = nullptr,
       int n_samples = 0,
       int n_features = 0,
       int leaf_size = 30,
       int dist_type = 1);
  Ball(const VectorVectorDouble &data, int leaf_size = 10, int dist_type = 1);
  Ball(const Db* db, int leaf_size = 30, int dist_type = 1, bool useSel = false);
  Ball(const Ball& p) = delete;
  Ball & operator=(const Ball& p) = delete;
  virtual ~Ball();

  KNN query(const double **test, int n_samples, int n_features, int n_neighbors = 1);
  KNN queryAsVVD(const VectorVectorDouble& test, int n_neighbors = 1);
  KNN queryOne(const double *test, int n_features, int n_neighbors = 1);
  KNN queryOneAsVD(const VectorDouble& test, int n_neighbors = 1);
  int queryClosest(const VectorDouble& test);
  void display(int level = -1) const;

protected:
  int _getFeatureNumber() const { return _tree->n_features; }
  int _getLeafSize() const { return _tree->leaf_size; }
  int _getSampleNumber() const { return _tree->n_samples; }
  int _getDistType() const { return _tree->dist_type; }

private:
  t_btree* _tree;
};
