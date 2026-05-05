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
#include "Tree/Ball.hpp"
#include "Tree/KNN.hpp"
#include "Tree/ball_algorithm.h"

namespace gstlrn
{
  class MeshEFaulted;

  class GSTLEARN_EXPORT BallFaulted: public Ball
  {
  public:
    BallFaulted() = default;

    BallFaulted(
      const MeshEFaulted* mesh,
      Id leaf_size = 10,
      bool all_available = true,
      Id default_distance_function = 1);

    ~BallFaulted() = default;

    VectorInt getIndices(const SpacePoint& Pt, Id n_neighbors = 1) override;
    Id queryOneInPlace(
      const VectorDouble& test,
      Id n_neighbors,
      VectorInt& indices,
      VectorDouble& distances,
      Id rank = 0) override;

  private:
    MatrixT<double> _getInformationFromMesh(
      const MeshEFaulted* mesh,
      Id* n_samples,
      Id* n_features);

    void _convertIndices(VectorInt& indices);

  private:
    VectorInt _mapIds;
  };
} // namespace gstlrn
