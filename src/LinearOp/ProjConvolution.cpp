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
#include "LinearOp/ProjConvolution.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorHelper.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "LinearOp/ProjMatrix.hpp"

#include "Matrix/MatrixSparse.hpp"
#include "Matrix/NF_Triplet.hpp"

ProjConvolution::ProjConvolution(const VectorDouble &convolution,
                                 const DbGrid *grid_point,
                                 const VectorInt& nodeRes2D,
                                 const VectorDouble& gext)
    : _convolution(convolution),
      _gridSeismic(grid_point),
      _nodeRes2D(nodeRes2D),
      _gext(gext),
      _shiftVector(),
      _gridSeis2D(nullptr),
      _gridRes2D(nullptr),
      _AProjHoriz(nullptr),
      _work()
{
  int ndim = grid_point->getNDim();
  if (ndim != 2 && ndim != 3)
  {
    messerr("ProjConvolution is limited to 2-D or 3-D case");
    return;
  }
  if (grid_point->getGrid().isRotated())
  {
    messerr("ProjConvolution is not implemented for Rotated grids yet");
    return;
  }

  _buildGridSeis2D();

  if (_nodeRes2D.empty())
    _nodeRes2D = _gridSeis2D->getNXs();

  _buildGridRes2D();

  _work.resize(_gridRes2D->getSampleNumber() * _gridSeismic->getNX(ndim-1));

  _buildAprojHoriz();

  _buildShiftVector();
}

ProjConvolution::~ProjConvolution()
{
  delete _gridSeis2D;
  delete _gridRes2D;
  if (_AProjHoriz != nullptr) delete _AProjHoriz;
}

void ProjConvolution::_buildGridSeis2D()
{
  int ndim = _getNDim();
  VectorInt nx_seis = _gridSeismic->getNXs();
  nx_seis.resize(ndim-1);
  VectorDouble dx_seis = _gridSeismic->getDXs();
  dx_seis.resize(ndim-1);
  VectorDouble x0_seis = _gridSeismic->getX0s();
  x0_seis.resize(ndim-1);
  _gridSeis2D = DbGrid::create(nx_seis,dx_seis,x0_seis);
}

void ProjConvolution::_buildGridRes2D()
{
  _gridRes2D = DbGrid::createCoveringDb(_gridSeis2D, _nodeRes2D, VectorDouble(),
                                        VectorDouble(), _gext);
}

int ProjConvolution::_buildAprojHoriz()
{
  // Create the Turbo Meshing on the Resolution 'ndim-1' grid
  MeshETurbo* mesh = MeshETurbo::createFromGrid(_gridRes2D);

  ProjMatrix* proj = ProjMatrix::create(_gridSeis2D, mesh);

  _AProjHoriz = proj->getAproj()->clone();

  delete mesh;
  delete proj;

  return 0;
}

/**
 * Calculate the vector of grid index shifts (in Point Grid)
 * This vector is calculated for the cell located in the center of the grid
 */
void ProjConvolution::_buildShiftVector()
{
  // Creating the characteristics of the Point Grid

  Grid grid = _getGridCharacteristicsRR();

  int ndim = _gridSeismic->getNDim();
  int center = 1;
  for (int idim = 0; idim < ndim; idim++)
    center *= grid.getNX(idim);
  center /= 2;

  VectorInt indp(ndim);
  VectorInt indm(ndim);
  _shiftVector.resize(_getConvSize());

  grid.rankToIndice(center, indp);
  for (int idim = 0; idim < ndim; idim++) indm[idim] = indp[idim];

  // Shift the index of last coordinate by the shift of the grid
  indp[ndim - 1] += _getHalfSize();

  for (int i = -_getHalfSize(); i <= _getHalfSize(); i++)
  {
    indm[ndim - 1] = indp[ndim - 1] + i;
    int id = grid.indiceToRank(indm);
    _shiftVector[i + _getHalfSize()] = id - center;
  }
}

bool ProjConvolution::_isVecDimCorrect(const VectorDouble &valonseismic,
                                       const VectorDouble &valonvertex) const
{
  if ((int) valonvertex.size() != getApexNumber())
  {
    messerr("Dimension of 'valonvertex'(%d) incorrect. If should be %d",
            (int) valonvertex.size(), getApexNumber());
    return false;
  }
  if ((int) valonseismic.size() != getPointNumber())
  {
    messerr("Dimension of 'valonseismic'(%d) incorrect. If should be %d",
            (int) valonseismic.size(), getPointNumber());
    return false;
  }
  if (_shiftVector.empty())
  {
    messerr("The ProjConvolution object has not been built correctly");
    return false;
  }
  return true;
}

/**
 * Apply the projection for a Seismic Grid Vector
 * and store the result on a Coarse Grid vector
 * @param valonseismic Input vector defined on the Seismic Grid
 * @param valonvertex  Output vector defined on the Coarse Grid
 * @return
 */
int ProjConvolution::point2mesh(const VectorDouble &valonseismic,
                                VectorDouble &valonvertex) const
{
  if (! _isVecDimCorrect(valonseismic, valonvertex)) return 1;

   int ndim  = _getNDim();

   // Get the characteristics of the R-R grid
   int slice_R = _gridRes2D->getSampleNumber();
   VectorDouble vec_R(slice_R);

   // Get the characteristics of the S-S grid
   int slice_S = _gridSeis2D->getSampleNumber();
   VectorDouble vec_S(slice_S);

   // Mesh barycenter on 'ndim-1' slices
   for (int iz = 0; iz < _gridSeismic->getNX(ndim-1); iz++)
   {
     VH::extractInPlace(valonseismic, vec_S, iz * slice_S);
     _AProjHoriz->prodMatVecInPlace(vec_S, vec_R, true);
     VH::mergeInPlace(vec_R, _work, iz * slice_R);
   }

   _convolveT(_work,valonvertex);
   return 0;
}

/**
 * Apply the Projection for a Coarse Grid vector
 * and store the result in a Seismic Grid Vector
 * @param valonvertex   Input vector defined on the Coarse Grid
 * @param valonseismic  Output vector defined on the Seismic grid
 * @return
 */
int ProjConvolution::mesh2point(const VectorDouble &valonvertex,
                                VectorDouble &valonseismic) const
{
  if (! _isVecDimCorrect(valonseismic, valonvertex)) return 1;

  int ndim  = _getNDim();

  // Get the characteristics of the R-R grid
  int slice_R = _gridRes2D->getSampleNumber();
  VectorDouble vec_R(slice_R);

  // Get the characteristics of the R-S grid
  int slice_S = _gridSeis2D->getSampleNumber();
  VectorDouble vec_S(slice_S);

  // Convolution
  _convolve(valonvertex, _work);

  // Mesh barycenter on 'ndim-1' slices
  for (int iz = 0; iz < _gridSeismic->getNX(ndim-1); iz++)
  {
    VH::extractInPlace(_work, vec_R, iz * slice_R);
    _AProjHoriz->prodMatVecInPlace(vec_R, vec_S, false);
    VH::mergeInPlace(vec_S, valonseismic, iz * slice_S);
  }

  return 0;
}

void ProjConvolution::_convolve(const VectorDouble &valonvertex,
                                VectorDouble &valonseismic) const
{

  for (auto &e : valonseismic)
     e = 0.;

  int count = (int) valonseismic.size();
  int size  = _getConvSize();
  double valp  = 0.;
  double valm = 0.;
  int id = 0;
  for (int is = 0; is < count; is++)
  {
    valp = 0;
    for (int j = 0; j < size; j++)
    {
      id = is + _shiftVector[j];
      valm = valonvertex[id];
      if( FFFF(valm))
      {
        valp = TEST;
        break;
      }
      valp += valm * _convolution[j];
    }
    valonseismic[is] = valp;
  }
}

void ProjConvolution::_convolveT(const VectorDouble &valonseismic,
                                 VectorDouble &valonvertex) const
{
  for (auto &e : valonvertex)
     e = 0.;

  int count = (int) valonseismic.size();
  int size  = _getConvSize();
  double valm = 0.;
  int id = 0;
  for (int is = 0; is < count; is++)
  {
    for (int j = 0; j < size; j++)
    {
      id = is + _shiftVector[j];
      valm = valonseismic[is];
      if (FFFF(valm))
      {
        valonvertex[id] = TEST;
        break;
      }
      valonvertex[id] += valm * _convolution[j];
    }
  }
}

Grid ProjConvolution::_getGridCharacteristicsRR(bool delLastDim) const
{
  int ndim = _gridSeismic->getNDim();

  VectorInt    nx = _gridRes2D->getNXs();
  VectorDouble dx = _gridRes2D->getDXs();
  VectorDouble x0 = _gridRes2D->getX0s();

  if (! delLastDim)
  {
    nx.resize(ndim);
    x0.resize(ndim);
    dx.resize(ndim);
    dx[ndim - 1] = _gridSeismic->getDX(ndim - 1);
    nx[ndim - 1] = _gridSeismic->getNX(ndim - 1) + (_getConvSize()-1);
    x0[ndim - 1] = _gridSeismic->getX0(ndim - 1) - (_getConvSize()-1) * dx[ndim-1];
  }
  Grid grid(ndim, nx, x0, dx);

  return grid;
}

/**
 * Grid matching Resolution in 'ndim-1' and Seismic for 'ndim'
 * @return
 */
Grid ProjConvolution::_getGridCharacteristicsRS() const
{
  int ndim = _gridSeismic->getNDim();

  Grid gridRR = _getGridCharacteristicsRR();
  VectorInt nxs   = gridRR.getNXs();
  VectorDouble dx = gridRR.getDXs();
  VectorDouble x0 = gridRR.getX0s();

  // Correct the last dimension
  nxs[ndim - 1] = _gridSeismic->getNX(ndim - 1);
  x0[ndim - 1]  = _gridSeismic->getX0(ndim - 1);
  Grid grid(ndim, nxs, x0, dx);

  return grid;
}

DbGrid* ProjConvolution::getResolutionGrid() const
{
  // Get the characteristics of the Point Grid
  Grid grid = _getGridCharacteristicsRR();

  // Create the new Point grid
  DbGrid* dbgrid = DbGrid::create(grid.getNXs(),
                                  grid.getDXs(),
                                  grid.getX0s());
  return dbgrid;
}

int ProjConvolution::getApexNumber() const
{
  Grid grid = _getGridCharacteristicsRR();
  return VH::product(grid.getNXs());
}

int ProjConvolution::getPointNumber() const
{
  VectorInt nxs = _gridSeismic->getNXs();
  return VH::product(nxs);
}
