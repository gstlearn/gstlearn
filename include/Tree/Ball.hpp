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
#include "Tree/ball.h"

class Db;

class GSTLEARN_EXPORT Ball
{
public:
  Ball(double **data = nullptr,
       int n_samples = 0,
       int n_features = 0,
       int leaf_size = 10,
       int dist_type = 0);
  Ball(VectorVectorDouble &data, int leaf_size = 10, int dist_type = 0);
  Ball(const Db* db, int leaf_size = 10, int dist_type = 0, bool useSel = false);
  Ball(const Ball& p);
  Ball & operator=(const Ball& p);
  virtual ~Ball();

  int getFeatureNumber() const { return _featureNumber; }
  int getLeafSize() const { return _leafSize; }
  int getSampleNumber() const { return _sampleNumber; }
  int getDistType() const { return _distType; }

  void setLeafSize(int leafSize) { _leafSize = leafSize; }
  void setData(double **data, int n_samples, int n_features);
  void setDataAsVVD(VectorVectorDouble& data);
  void setDistType(int distType) { _distType = distType; }

  int build() const;
  t_knn query(double **test, int n_samples, int n_features, int n_neighbors);
  t_knn queryAsVVD(VectorVectorDouble& test, int n_neighbors);
  t_knn queryOne(double *test, int n_features, int n_neighbors);
  t_knn queryOne(VectorDouble& test, int n_neighbors);

private:
  int _leafSize;
  int _sampleNumber;
  int _featureNumber;
  int _distType;
  double **_data;

  mutable t_btree* _tree;
  mutable double (* _distFunction)(double, double , int);
};
