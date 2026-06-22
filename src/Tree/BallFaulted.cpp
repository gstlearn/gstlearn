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

#include "Tree/BallFaulted.hpp"
#include "Basic/VectorHelper.hpp"
#include "Mesh/MeshEFaulted.hpp"
#include "Space/SpacePoint.hpp"
#include "Tree/Ball.hpp"
#include "Tree/ball_algorithm.h"

#include "geoslib_define.h"

namespace gstlrn
{
  /**
   * @brief Construct a new BallFaulted tree object based on the barycenters of the meshes
   *
   * @param mesh  AMesh description
   * @param leaf_size Number of elements in the leafs of the BallFaulted tree
   * @param all_available True if samples are available for selection at the beginning
   * @param default_distance_function 1 for Euclidean distance, 2 for Manhattan
   */
  BallFaulted::BallFaulted(
    const MeshEFaulted* mesh,
    Id leaf_size,
    bool all_available,
    Id default_distance_function)
    : Ball()
    , _mapIds()
  {
    Id n_samples;
    Id n_features;
    auto internal = _getInformationFromMesh(mesh, &n_samples, &n_features);
    if (internal.empty()) return;

    _tree = t_btree(
      std::move(internal), n_samples, n_features, all_available, leaf_size,
      default_distance_function);
  }

  VectorInt BallFaulted::getIndices(const SpacePoint& Pt, Id n_neighbors)
  {
    VectorInt indices = Ball::getIndices(Pt, n_neighbors);
    _convertIndices(indices);

    return indices;
  }

  Id BallFaulted::queryOneInPlace(
    const VectorDouble& test,
    Id n_neighbors,
    VectorInt& indices,
    VectorDouble& distances,
    Id rank)
  {

    Id out = Ball::queryOneInPlace(test, n_neighbors, indices, distances, rank);
    _convertIndices(indices);
    return out;
  }

  MatrixT<double> BallFaulted::_getInformationFromMesh(
    const MeshEFaulted* mesh,
    Id* n_samples,
    Id* n_features)
  {
    VectorDouble oneColumn;
    Id ndim = mesh->getNDim();
    Id nmesh = mesh->getNMeshes();
    Id nFaultMeshes = mesh->getNVirtualFaultMeshes();
    VectorInt faultIds = mesh->getFaultIds();

    // Core allocation
    MatrixT<double> internal(nmesh - nFaultMeshes, ndim);
    _mapIds.reserve(nmesh - nFaultMeshes);

    Id shift = 0;
    // Loading the information from mesh
    for (Id imesh = 0; imesh < nmesh; imesh++)
    {
      if (faultIds[imesh] == 0)
      {
        mesh->getBarycenterInPlace(imesh, internal.getRow(imesh - shift));
        _mapIds.insert(imesh - shift, imesh);
      }
      else
        ++shift;
    }

    *n_samples = nmesh - nFaultMeshes;
    *n_features = ndim;
    return internal;
  }

  void BallFaulted::_convertIndices(VectorInt& indices)
  {
    for (size_t i_indice = 0; i_indice < indices.size(); ++i_indice)
    {
      indices[i_indice] = _mapIds[indices[i_indice]];
    }
  }
} // namespace gstlrn
