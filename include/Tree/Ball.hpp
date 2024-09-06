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
class SpacePoint;

class GSTLEARN_EXPORT Ball
{
public:
  Ball(const double** data               = nullptr,
       int n_samples                     = 0,
       int n_features                    = 0,
       int leaf_size                     = 10,
       double (*dist_function)(const double* x1,
                               const double* x2,
                               int size) = nullptr);
  Ball(const VectorVectorDouble& data,
       double (*dist_function)(const double* x1,
                               const double* x2,
                               int size) = nullptr,
       int leaf_size                     = 10);
  Ball(const Db* db,
       double (*dist_function)(const double* x1,
                               const double* x2,
                               int size) = nullptr,
       int leaf_size                     = 10,
       bool useSel                       = false);
  Ball(const Ball& p)            = delete;
  Ball& operator=(const Ball& p) = delete;
  virtual ~Ball();

  void init(const Db* db,
            double (*dist_function)(const double* x1,
                                    const double* x2,
                                    int size) = nullptr,
            int leaf_size                     = 10,
            bool useSel                       = false);

  KNN query(const double** test,
            int n_samples,
            int n_features,
            int n_neighbors = 1);
  KNN queryAsVVD(const VectorVectorDouble& test, int n_neighbors = 1);
  KNN queryOne(const double *test, int n_features, int n_neighbors = 1);
  KNN queryOneAsVD(const VectorDouble& test, int n_neighbors = 1);
  KNN queryOneAsVD(const SpacePoint& Pt, int n_neighbors = 1);
  VectorInt getIndices(const SpacePoint& Pt, int n_neighbors = 1);
  int queryClosest(const VectorDouble& test);
  void display(int level = -1) const;

protected:
  int _getFeatureNumber() const { return _tree->n_features; }
  int _getLeafSize() const { return _tree->leaf_size; }
  int _getSampleNumber() const { return _tree->n_samples; }

private:
  t_btree* _tree;
};
