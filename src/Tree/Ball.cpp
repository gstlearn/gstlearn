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
#include "Tree/ball.h"
#include "Db/Db.hpp"
#include "Basic/AStringable.hpp"

Ball::Ball(double **data, int n_samples, int n_features, int leaf_size, int dist_type)
    : _leafSize(leaf_size),
      _sampleNumber(n_samples),
      _featureNumber(n_features),
      _distType(dist_type),
      _data(nullptr),
      _tree(nullptr),
      _distFunction()
{
  _data = copy_double_arr(data, _sampleNumber, _featureNumber);
}

Ball::Ball(VectorVectorDouble& data, int leaf_size, int dist_type)
    : _leafSize(leaf_size),
      _sampleNumber(0),
      _featureNumber(0),
      _distType(dist_type),
      _data(nullptr),
      _tree(nullptr),
      _distFunction()
{
  _sampleNumber = (int) data[0].size();
  _featureNumber = (int) data.size();
  _data = copy_double_arrAsVVD(data);
}

Ball::Ball(const Db* db, int leaf_size, int dist_type, bool useSel)
: _leafSize(leaf_size),
  _sampleNumber(0),
  _featureNumber(0),
  _distType(dist_type),
  _data(nullptr),
  _tree(nullptr),
  _distFunction()
{
  VectorVectorDouble data = db->getAllCoordinates(useSel);
  _sampleNumber = (int) data[0].size();
  _featureNumber = (int) data.size();
  _data = copy_double_arrAsVVD(data);
}

Ball::Ball(const Ball &p)
    : _leafSize(p._leafSize),
      _sampleNumber(p._sampleNumber),
      _featureNumber(p._featureNumber),
      _distType(p._distType),
      _data(nullptr),
      _tree(nullptr),_distFunction()
{
  _data = copy_double_arr(p._data, p._sampleNumber, p._featureNumber);

  // TODO: the _tree is not copied (function does not exist)
  _tree = nullptr;
}

Ball& Ball::operator=(const Ball& p)
{
  if (this != &p)
  {
    _leafSize = p._leafSize;
    _sampleNumber = p._sampleNumber;
    _featureNumber = p._featureNumber;
    _distType = p._distType;

    _data = copy_double_arr(p._data, p._sampleNumber, p._featureNumber);

    // TODO: the _tree is not copied (function does not exist)
    _tree = nullptr;

    _distFunction = p._distFunction;
  }
  return *this;
}

Ball::~Ball()
{
  free_2d_double(_data, _sampleNumber);
  free_tree(_tree);
}

void Ball::setData(double **data, int n_samples, int n_features)
{
  _sampleNumber = n_samples;
  _featureNumber = n_features;
  _data = copy_double_arr(data, _sampleNumber, _featureNumber);
}

void Ball::setDataAsVVD(VectorVectorDouble& data)
{
  _sampleNumber = (int) data[0].size();
  _featureNumber = (int) data.size();
  _data = copy_double_arrAsVVD(data);
}

int Ball::build() const
{
  if (_leafSize <= 0)
  {
    messerr("You must define 'leafSize' beforehand");
    return 1;
  }
  if (_data == nullptr) return 1;

  // Clean already existing tree (if any)

  free_tree(_tree);

  _tree = btree_init(_data, _sampleNumber, _featureNumber, _leafSize, _distType);
  return 0;
}

t_knn Ball::query(double **test, int n_samples, int n_features, int n_neighbors)
{
  t_knn knn;
  if (_tree == nullptr)
  {
    messerr("Tree must be 'built' beforehand");
    return knn;
  }
  if (test == nullptr) return knn;
  if (n_features != _featureNumber)
  {
    messerr("Tree (%d) and Test (%d) must share the same feature dimension",
            n_features, _featureNumber);
    return knn;
  }
  knn = btree_query(_tree, test, n_samples, n_features, n_neighbors);
  return knn;
}

t_knn Ball::queryAsVVD(VectorVectorDouble& test, int n_neighbors)
{
  t_knn knn;
  if (_tree == nullptr)
  {
    messerr("Tree must be 'built' beforehand");
    return knn;
  }
  if (test.empty()) return knn;
  int n_samples = (int) test[0].size();
  int n_features = (int) test.size();
  if (n_features != _featureNumber)
  {
    messerr("Tree (%d) and Test (%d) must share the same feature dimension",
            n_features, _featureNumber);
    return knn;
  }
  knn = btree_query(_tree, copy_double_arrAsVVD(test), n_samples, n_features, n_neighbors);
  return knn;
}

t_knn Ball::queryOne(double *test, int n_features, int n_neighbors)
{
  t_knn knn;
  if (_tree == nullptr)
  {
    messerr("Tree must be 'built' beforehand");
    return knn;
  }
  if (test == nullptr) return knn;
  if (n_features != _featureNumber)
  {
    messerr("Tree (%d) and Test (%d) must share the same feature dimension",
            n_features, _featureNumber);
    return knn;
  }
  double** internal = &test;
  knn = btree_query(_tree, internal, 1, n_features, n_neighbors);
  return knn;
}

t_knn Ball::queryOneAsVD(VectorDouble& test, int n_neighbors)
{
  t_knn knn;
  if (_tree == nullptr)
  {
    messerr("Tree must be 'built' beforehand");
    return knn;
  }
  if (test.empty()) return knn;
  int n_features = (int) test.size();
  if (n_features != _featureNumber)
  {
    messerr("Tree (%d) and Test (%d) must share the same feature dimension",
            n_features, _featureNumber);
    return knn;
  }
  double* local = test.data();
  double **internal = &local;
  knn = btree_query(_tree, internal, 1, n_features, n_neighbors);
  return knn;
}


void display(t_knn& knn, int ns_max, int nn_max)
{
  int ns = knn.n_samples;
  if (ns_max >= 0) ns = MIN(ns, ns_max);
  int nn = knn.n_neighbors;
  if (nn_max > 0) nn = MIN(nn, nn_max);
  for (int is = 0; is < ns; is++)
  {
    message("Indices = ");
    for (int in = 0; in < nn; in++)
      message(" %d", knn.indices[is][in]);
    message("\n");

    message("Distances = ");
    for (int in = 0; in < nn; in++)
      message(" %lf", knn.distances[is][in]);
    message("\n");
  }
}

VectorInt getIndices(t_knn& knn, int rank)
{
  if (rank < 0 || rank >= knn.n_samples) return VectorInt();
  return VectorHelper::initVInt(knn.indices[rank], knn.n_neighbors);
}
VectorDouble getDistance(t_knn& knn, int rank)
{
  if (rank < 0 || rank >= knn.n_samples) return VectorDouble();
  return VectorHelper::initVDouble(knn.distances[rank], knn.n_neighbors);
}
