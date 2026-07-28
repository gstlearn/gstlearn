/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "LinearOp/ShiftOpStencil.hpp"
#include "Basic/Grid.hpp"
#include "Basic/Indirection.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "LinearOp/AShiftOp.hpp"
#include "LinearOp/ShiftOpMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "geoslib_define.h"

#include <omp.h>

namespace gstlrn
{
  ShiftOpStencil::ShiftOpStencil(
    const MeshETurbo* mesh,
    const CovAniso* cova,
    bool verbose)
    : AShiftOp()
    , _stencil()
    , _isInside()
    , _useLambdaSingleVal(true)
    , _mesh()
  {
    if (_buildInternal(mesh, cova, verbose)) return;
  }

  Id ShiftOpStencil::_addToDest(const constvect inv, vect outv) const
  {
    const VectorDouble& currentWeights = _stencil.weights();

    auto nw = _stencil.getNWeights();
    Id size = static_cast<Id>(inv.size());
    const Indirection& indirect = _mesh->getGridIndirect();

    auto nbthread = static_cast<I32>(OptCustom::query("ompthreads", 1));

    if (!indirect.isDefined())
    {
      const VectorInt& absoluteShifts = _stencil.absoluteShifts();
      // Use the fast option when no selection is defined on the Grid
#pragma omp parallel for num_threads(nbthread)
      for (Id ic = 0; ic < size; ic++)
      {
        double total = 0.;
        if (_isInside[ic])
        {
          for (Id iw = 0; iw < nw; iw++)
          {
            Id iabs = ic + absoluteShifts[iw];
            if (_isInside[iabs])
            {
              total += currentWeights[iw] * inv[iabs];
            }
          }
        }
        outv[ic] = total;
      }
    }
    else
    {
      const Grid& grid = _mesh->getGrid();
      Id ndim = _mesh->getNDim();
      const VectorVectorInt& relativeShifts = _stencil.relativeShifts();

      VectorInt center(ndim);
      VectorInt local(ndim);

#pragma omp parallel for firstprivate(center, local) num_threads(nbthread)
      for (Id ic = 0; ic < size; ic++)
      {
        double total = 0.;

        // Check if the target point is not on the edge and not masked
        if (_isInside[ic])
        {
          Id rank = indirect.getRToA(ic);
          grid.rankToIndice(rank, center);
          for (Id iw = 0; iw < nw; iw++)
          {
            VH::add(local, center, relativeShifts[iw]);
            Id ie = grid.indiceToRank(local);
            rank = indirect.getAToR(ie);
            if (rank >= 0) total += currentWeights[iw] * inv[rank];
          }
        }
        outv[ic] = total;
      }
    }
    return 0;
  }

  void ShiftOpStencil::resetModif() const
  {
    _stencil.resetModif();
  }

  void ShiftOpStencil::normalizeLambdaBySills(const AMesh* mesh)
  {
    if (_cova->isNoStatForVariance())
    {
      _Lambda.resize(_napices);
      for (auto& e: _Lambda)
      {
        e = _stencil.getLambda();
      }
      AShiftOp::normalizeLambdaBySills(mesh);
      _useLambdaSingleVal = false;
    }
    else
    {
      _stencil.getLambda() /= sqrt(_cova->getSill(0, 0));
    }
  }

  double ShiftOpStencil::_getMaxEigenValue() const
  {
    double s = 0.;
    for (const auto& e: _stencil.weights())
    {
      s += ABS(e);
    }
    return s;
  }

  double ShiftOpStencil::getLambda(Id iapex) const
  {
    if (_useLambdaSingleVal) return _stencil.getLambda();
    return AShiftOp::getLambda(iapex);
  }

  double ShiftOpStencil::logDetLambda() const
  {
    if (_useLambdaSingleVal) return 2. * log(_stencil.getLambda()) * _napices;
    return AShiftOp::logDetLambda();
  }

  void ShiftOpStencil::multiplyByValueAndAddDiagonal(double v1, double v2) const
  {
    _stencil.multiplyByValueAndAddDiagonal(v1, v2);
  }

  Id ShiftOpStencil::_buildInternal(
    const MeshETurbo* mesh,
    const CovAniso* cova,
    bool verbose)
  {
    // Preliminary checks

    if (cova == nullptr)
    {
      messerr("The argument 'cova' must be provided");
      return 1;
    }
    if (mesh == nullptr)
    {
      messerr("The argument 'mesh' must be provided");
      return 1;
    }

    _mesh = mesh;
    _napices = mesh->getNApices();
    Id ndim = _mesh->getNDim();

    _setCovAniso(cova);

    if (_cova->isNoStatForAnisotropy())
    {
      messerr("The Shiftop as a Stencil is incompatible with non-stationarity");
      return 1;
    }

    _stencil = ShiftStencil(_mesh, _cova.get(), verbose);

    const Grid& grid = mesh->getGrid();
    VectorInt center(ndim);
    VectorInt NXs = grid.getNXs();

    // Delineate the border of the grid (not to be treated)
    Id size = _mesh->getNApices();
    _isInside.fill(true, size);

    const Indirection& indirect = _mesh->getGridIndirect();
    for (Id i = 0; i < size; i++)
    {
      Id rank = indirect.isDefined() ? indirect.getRToA(i) : i;
      grid.rankToIndice(rank, center);
      bool flagInside = true;
      for (Id idim = 0; idim < ndim && flagInside; idim++)
      {
        Id ival = center[idim];
        if (ival <= 0 || ival >= NXs[idim] - 1) flagInside = false;
      }
      _isInside[i] = flagInside;
    }

    if (verbose)
    {
      Id ntreated = 0;
      for (Id i = 0; i < size; i++)
        if (_isInside[i]) ntreated++;
      message(
        "Number of pixels inside the grid (no edge effect) = %d/%d\n", ntreated,
        size);
    }

    return 0;
  }

  ShiftStencil::ShiftStencil(
    const MeshETurbo* mesh,
    const CovAniso* cova,
    bool verbose)
  {
    // Create a local Turbo Meshing starting from a DbGrid (dimension 3 should be sufficient but 5 is safer)
    Id ndim = mesh->getNDim();
    const Grid& grid = mesh->getGrid();
    VectorInt NXs = grid.getNXs();
    VectorInt nxlocal(ndim, 5);

    MeshETurbo localMesh(
      nxlocal, grid.getDXs(), grid.getX0s(), grid.getRotAngles());
    ShiftOpMatrix shiftMat(&localMesh, cova, nullptr, verbose);

    // Display the vector of the 'S' matrix for the center Apex
    MatrixSparse* S = shiftMat.getS();
    auto centerApex = localMesh.getNApices() / 2;
    VectorDouble centerColumn = S->getColumn(centerApex);

    // Fill lambda
    _lambdaVal = shiftMat.getLambda(centerApex);

    // Get the indices of the centerApex
    VectorInt center(ndim);
    VectorInt other(ndim);
    localMesh.getApexIndicesInPlace(centerApex, center);

    // Get the non-zero elements of the center column
    for (Id i = 0, n = static_cast<Id>(centerColumn.size()); i < n; i++)
    {
      double weight = centerColumn[i];
      if (ABS(weight) < EPSILON6) continue;
      localMesh.getApexIndicesInPlace(i, other);
      other -= center;
      _relativeShifts.push_back(other);
      _weights.push_back(weight);
    }

    auto nw = getNWeights();

    // Calculate the shifts (from the center cell) for each weight
    // This is calculated for a reference pixel (center of the grid)
    _absoluteShifts.fill(0, nw);
    for (Id idim = 0; idim < ndim; idim++) center[idim] = NXs[idim] / 2;
    Id iorigin = grid.indiceToRank(center);

    VectorInt local(ndim);
    for (Id iw = 0; iw < nw; iw++)
    {
      VH::add(local, center, _relativeShifts[iw]);
      Id iabs = grid.indiceToRank(local);
      _absoluteShifts[iw] = iabs - iorigin;
    }

    // Print the contents of non-zero elements
    if (verbose) print();
  }

  const VectorDouble& ShiftStencil::weights() const
  {
    if (_useModifiedShift)
    {
      return _weightsSimu;
    }
    else
    {
      return _weights;
    }
  }

  void ShiftStencil::multiplyByValueAndAddDiagonal(double v1, double v2) const
  {
    _weightsSimu = VectorDouble(_weights.size());
    for (Id i = 0; i < static_cast<Id>(_weights.size()); i++)
      _weightsSimu[i] = v1 * _weights[i];
    Id center = static_cast<Id>(_relativeShifts.size()) / 2;
    _weightsSimu[center] += v2;
    _useModifiedShift = true;
  }

  void ShiftStencil::resetModif() const
  {
    _useModifiedShift = false;
  }

  void ShiftStencil::print() const
  {
    auto nweight = _weights.size();
    Id ndim = _relativeShifts.front().size();

    mestitle(0, "Stencil contents");
    for (size_t i = 0; i < nweight; i++)
    {
      message("Weight %d/%d - Relative (", i + 1, nweight);
      for (Id idim = 0; idim < ndim; idim++)
        message("%2d ", _relativeShifts[i][idim]);
      message(") - Absolute (%4d)", _absoluteShifts[i]);
      message(" : %lf\n", _weights[i]);
    }
  }
} // namespace gstlrn
