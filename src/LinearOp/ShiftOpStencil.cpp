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
#include "LinearOp/ShiftOpMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Covariances/CovAniso.hpp"
#include "Basic/Grid.hpp"

#include "geoslib_define.h"

ShiftOpStencil::ShiftOpStencil(const MeshETurbo* mesh,
                               const CovAniso* cova,
                               bool verbose)
  : AShiftOp()
  , _ndim(0)
  , _relativeRanks()
  , _weights()
  , _locations()
  , _onEdge()
{
  _buildInternal(mesh, cova, verbose);
}

ShiftOpStencil::ShiftOpStencil(const ShiftOpStencil& shift)
  : AShiftOp(shift)
  , _ndim(shift._ndim)
  , _relativeRanks(shift._relativeRanks)
  , _weights(shift._weights)
  , _locations(shift._locations)
  , _onEdge(shift._onEdge)
{
}

ShiftOpStencil& ShiftOpStencil::operator=(const ShiftOpStencil& shift)
{
  if (this != &shift)
  {
    _ndim          = shift._ndim;
    _relativeRanks = shift._relativeRanks;
    _weights       = shift._weights;
    _locations     = shift._locations;
    _onEdge        = shift._onEdge;
  }
  return *this;
}

ShiftOpStencil::~ShiftOpStencil() {}

int ShiftOpStencil::_addToDest(const constvect inv, vect outv) const
{
  DECLARE_UNUSED(inv, outv)
  return 0;
}

void ShiftOpStencil::normalizeLambdaBySills(const AMesh* mesh)
{
  DECLARE_UNUSED(mesh)
}

void ShiftOpStencil::prodLambda(const constvect x,
                                vect y,
                                const EPowerPT& power) const
{
  DECLARE_UNUSED(x, y, power)
}

double ShiftOpStencil::getMaxEigenValue() const
{
  return TEST;
}

void ShiftOpStencil::_buildInternal(const MeshETurbo* mesh,
                                    const CovAniso* cova,
                                    bool verbose)
{
  _ndim = mesh->getNDim();

  // Create a local Turbo Meshing starting from a DbGrid

  const Grid& grid = mesh->getGrid();
  VectorInt nx     = grid.getNXs();
  for (int idim = 0; idim < _ndim; idim++) nx[idim] = 5;

  MeshETurbo localMesh(nx, grid.getDXs(), grid.getX0s(), grid.getRotAngles());
  localMesh.display();

  ShiftOpMatrix shiftMat(&localMesh, cova, nullptr, verbose);

  // Display the vector of the 'S' matrix for the center Apex
  MatrixSparse* S           = shiftMat.getS();
  int centerApex            = localMesh.getNApices() / 2;
  VectorDouble centerColumn = S->getColumn(centerApex);

  // Get the indices of the centerApex
  VectorInt center(_ndim);
  VectorInt other(_ndim);
  localMesh.getApexIndicesInPlace(centerApex, center);
  VH::display("Ranks of Center Apex", center);

  // Get the non-zero elements of the center column
  _relativeRanks.clear();
  _weights.clear();
  _locations.clear();
  for (int i = 0, n = (int)centerColumn.size(); i < n; i++)
  {
    double weight = centerColumn[i];
    if (ABS(weight) < EPSILON6) continue;
    localMesh.getApexIndicesInPlace(i, other);
    VH::subtractInPlace(other, center);
    _relativeRanks.push_back(other);
    _weights.push_back(weight);
  }

  // Print the contents of non-zero elements
  if (verbose) _printStencil();
}

void ShiftOpStencil::_printStencil() const
{
  int nweight = (int)_weights.size();
  mestitle(0, "Stencil contents");
  for (int i = 0; i < nweight; i++)
  {
    message("Weight %d/%d (", i + 1, nweight);
    for (int idim = 0; idim < _ndim; idim++)
      message("%d ", _relativeRanks[i][idim]);
    message(") : %lf\n", _weights[i]);
  }
}