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

#include "Basic/AStringable.hpp"
#include "Basic/Indirection.hpp"
#include "LinearOp/ShiftOpMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Covariances/CovAniso.hpp"
#include "Basic/Grid.hpp"

#include "geoslib_define.h"

ShiftOpStencil::ShiftOpStencil(const MeshETurbo* mesh,
                               const CovAniso* cova,
                               bool verbose)
  : AShiftOp()
  , _relativeShifts()
  , _absoluteShifts()
  , _weights()
  , _isInside()
  , _mesh()
{
  if (_buildInternal(mesh, cova, verbose)) return;
}

ShiftOpStencil::ShiftOpStencil(const ShiftOpStencil& shift)
  : AShiftOp(shift)
  , _relativeShifts(shift._relativeShifts)
  , _absoluteShifts(shift._absoluteShifts)
  , _weights(shift._weights)
  , _isInside(shift._isInside)
  , _mesh(shift._mesh)
{
}

ShiftOpStencil& ShiftOpStencil::operator=(const ShiftOpStencil& shift)
{
  if (this != &shift)
  {
    _relativeShifts = shift._relativeShifts;
    _absoluteShifts = shift._absoluteShifts;
    _weights        = shift._weights;
    _isInside       = shift._isInside;
    _mesh           = shift._mesh;
  }
  return *this;
}

ShiftOpStencil::~ShiftOpStencil() {}

int ShiftOpStencil::_addToDest(const constvect inv, vect outv) const
{
  int nw = _getNWeights();
  int size = (int) inv.size();
  const Indirection& indirect = _mesh->getGridIndirect();

  double total;
  if (!indirect.isDefined())
  {
    for (int ic = 0; ic < size; ic++)
    {
      total = 0.;
      if (_isInside[ic])
      {
        for (int iw = 0; iw < nw; iw++)
        {
          int iabs = ic + _absoluteShifts[iw];
          total += _weights[iw] * inv[iabs];
        }
      }
      outv[ic] = total;
    }
  }
  else
  {
    const Grid& grid = _mesh->getGrid();
    int ndim         = _mesh->getNDim();

    VectorInt center(ndim);
    VectorInt local(ndim);
    for (int ic = 0; ic < size; ic++)
    {
      total = 0.;

      // Check if the target point is not on the edge and not masked
      if (_isInside[ic] && indirect.getAToR(ic) >= 0)
      {
        grid.rankToIndice(ic, center);
        for (int iw = 0; iw < nw; iw++)
        {
          local = center;
          VH::addInPlace(local, _relativeShifts[iw]);
          int ie = grid.indiceToRank(local);
          if (indirect.getAToR(ie) >= 0) total += _weights[iw] * inv[ie];
        }
      }
      outv[ic] = total;
    }
    }
    return 0;
  }

void ShiftOpStencil::normalizeLambdaBySills(const AMesh* mesh)
{
  DECLARE_UNUSED(mesh)
}

void ShiftOpStencil::addProdLambda(const constvect x,
                                   vect y,
                                   const EPowerPT& power) const
{
  DECLARE_UNUSED(x, y, power)
}

double ShiftOpStencil::getMaxEigenValue() const
{
  return TEST;
}

int ShiftOpStencil::_buildInternal(const MeshETurbo* mesh,
                                   const CovAniso* cova,
                                   bool verbose)
{
  _mesh = mesh;
  int ndim        = _mesh->getNDim();

  // Preliminary checks
  if (cova->isNoStat())
  {
    messerr("The Shiftp as a Stencil is incompatible with non-stationarity");
    return 1;
  }

  // Create a local Turbo Meshing starting from a DbGrid (Dimension 5 sgould be sufficient)
  const Grid& grid = mesh->getGrid();
  VectorInt NXs    = grid.getNXs();
  VectorInt nxlocal(ndim, 5);

  MeshETurbo localMesh(nxlocal, grid.getDXs(), grid.getX0s(), grid.getRotAngles());
  ShiftOpMatrix shiftMat(&localMesh, cova, nullptr, verbose);

  // Display the vector of the 'S' matrix for the center Apex
  MatrixSparse* S           = shiftMat.getS();
  int centerApex            = localMesh.getNApices() / 2;
  VectorDouble centerColumn = S->getColumn(centerApex);

  // Get the indices of the centerApex
  VectorInt center(ndim);
  VectorInt other(ndim);
  localMesh.getApexIndicesInPlace(centerApex, center);

  // Get the non-zero elements of the center column
  _relativeShifts.clear();
  _weights.clear();
  for (int i = 0, n = (int)centerColumn.size(); i < n; i++)
  {
    double weight = centerColumn[i];
    if (ABS(weight) < EPSILON6) continue;
    localMesh.getApexIndicesInPlace(i, other);
    VH::subtractInPlace(other, center);
    _relativeShifts.push_back(other);
    _weights.push_back(weight);
  }
  int nw = _getNWeights();

  // Calculate the shifts (from the center cell) for each weight
  // This is calculated for a reference pixel (center of the grid)
  _absoluteShifts.fill(0, nw);
  for (int idim = 0; idim < ndim; idim++) center[idim] = NXs[idim] / 2;
  int iorigin = grid.indiceToRank(center);

  for (int iw = 0; iw < nw; iw++)
  {
    VectorInt local = center;
    VH::addInPlace(local, _relativeShifts[iw]);
    int iabs = grid.indiceToRank(local);
    _absoluteShifts[iw] = iabs - iorigin;
  }

  // Delineate the border of the grid (not to be treated)
  int size = _mesh->getNApices();
  _isInside.fill(true, size);

  for (int i = 0; i < size; i++)
  {
    grid.rankToIndice(i, center);
    bool flagInside = true;
    for (int idim = 0; idim < ndim && flagInside; idim++)
    {
      int ival = center[idim];
      if (ival <= 0 || ival >= NXs[idim] - 1) flagInside = false;
    }
    _isInside[i] = flagInside;
  }

  // Print the contents of non-zero elements
  if (verbose) _printStencil();

  return 0;
}

void ShiftOpStencil::_printStencil() const
{
  int nweight = _getNWeights();
  int ndim    = _mesh->getNDim();
  int size    = _mesh->getNApices();
  mestitle(0, "Stencil contents");
  for (int i = 0; i < nweight; i++)
  {
    message("Weight %d/%d - Relative (", i + 1, nweight);
    for (int idim = 0; idim < ndim; idim++)
      message("%2d ", _relativeShifts[i][idim]);
    message(") - Absolute (%4d)", _absoluteShifts[i]);
    message(" : %lf\n", _weights[i]);
  }

  int ntreated = 0;
  for (int i = 0; i < size; i++)
    if (_isInside[i]) ntreated++;
  message("Number of pixels inside the grid (no edge effect) = %d/%d\n", ntreated, size);
}