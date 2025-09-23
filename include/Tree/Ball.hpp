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

#include "Matrix/MatrixT.hpp"
#include "Tree/KNN.hpp"
#include "Tree/ball_algorithm.h"

namespace gstlrn
{
class Db;
class AMesh;
class SpacePoint;

class GSTLEARN_EXPORT Ball
{
public:
  Ball() = default;

  Ball(const Db* dbin,
       const Db* dbout              = nullptr,
       Id leaf_size                 = 10,
       bool all_available           = true,
       Id default_distance_function = 1,
       bool useSel                  = false);
  Ball(const AMesh* mesh,
       Id leaf_size                 = 10,
       bool all_available           = true,
       Id default_distance_function = 1);

  void init(const Db* db,
            Id leaf_size                 = 10,
            bool all_available           = true,
            Id default_distance_function = 1,
            bool useSel                  = false);

  KNN queryAsVVD(const VectorVectorDouble& test, Id n_neighbors = 1);
  KNN queryOne(const double* test, Id n_features, Id n_neighbors = 1);
  KNN queryOneAsVD(const VectorDouble& test, Id n_neighbors = 1);
  KNN queryOneAsVDFromSP(const SpacePoint& Pt, Id n_neighbors = 1);
  VectorInt getIndices(const SpacePoint& Pt, Id n_neighbors = 1);
  Id queryClosest(const VectorDouble& test);
  Id queryOneInPlace(const VectorDouble& test,
                     Id n_neighbors,
                     VectorInt& indices,
                     VectorDouble& distances,
                     Id rank = 0);
  void display(Id level = -1) const;
  Id setAvailable(Id rank, bool status);
  Id resetAvailable(bool status);
  bool empty() const { return _tree.data.empty(); }

protected:
  Id _getFeatureNumber() const { return _tree.n_features; }
  Id _getLeafSize() const { return _tree.leaf_size; }
  Id _getNSample() const { return _tree.n_samples; }

private:
  bool _isAvailableDefined() const;
  static MatrixT<double> _getInformationFromDb(const Db* dbin,
                                               const Db* dbout,
                                               bool useSel,
                                               Id* n_samples,
                                               Id* n_features);
  static MatrixT<double> _getInformationFromMesh(const AMesh* mesh,
                                                 Id* n_samples,
                                                 Id* n_features);

private:
  t_btree _tree;
};

GSTLEARN_EXPORT MatrixT<Id> findNN(const Db* dbin,
                                   const Db* dbout              = nullptr,
                                   Id nb_neigh                  = 3,
                                   bool flagShuffle             = false,
                                   bool verbose                 = false,
                                   Id leaf_size                 = 10,
                                   Id default_distance_function = 1);
} // namespace gstlrn
