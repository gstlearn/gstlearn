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
#include "Basic/Grid.hpp"

#include "Geometry/Rotation.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "geoslib_define.h"

#include <math.h>

typedef struct
{
  int curech;
  int ndim;
  VectorInt nx;
  VectorInt order;
  VectorInt indg;
  VectorInt tab;
} DimLoop;

/****************************************************************************/
/*!
 ** Perform the recursion through the dimension
 **
 ** \param[in]  idim    Space rank
 ** \param[in]  verbose Verbose flag
 ** \param[in]  dlp     DimLoop structure
 **
 *****************************************************************************/
static void _dimensionRecursion(int idim, bool verbose, DimLoop& dlp)
{
  int order, sdim, nval, ival, ndim;

  // Assignments

  ndim = dlp.ndim;

  if (idim < 0)
  {

    // We have reached the bottom of the pile, evaluate the absolute address

    ival = dlp.indg[ndim - 1];
    for (int jdim = ndim - 2; jdim >= 0; jdim--)
      ival = ival * dlp.nx[jdim] + dlp.indg[jdim];
    dlp.tab[dlp.curech++] = ival + 1;

    // Optional printout

    if (verbose)
    {
      message("node (");
      for (int jdim = 0; jdim < ndim; jdim++)
        message(" %d", dlp.indg[jdim] + 1);
      message(" ) -> %d\n", ival + 1);
    }
    return;
  }

  order = dlp.order[idim];
  sdim = ABS(order) - 1;
  nval = dlp.nx[sdim];

  // Loop

  for (int jy = 0; jy < nval; jy++)
  {
    dlp.indg[sdim] = (order < 0) ? nval - jy - 1 : jy;
    _dimensionRecursion(idim - 1, verbose, dlp);
  }
}

Grid::Grid(int ndim,
           const VectorInt &nx,
           const VectorDouble &x0,
           const VectorDouble &dx)
  : AStringable(),
    _nDim(ndim)
  , _nx()
  , _x0()
  , _dx()
  , _rotation()
  , _iter(0)
  , _nprod(0)
  , _counts()
  , _order()
  , _indices()
  , _iwork0(ndim)
  , _work1(ndim)
  , _work2(ndim)
{
  _allocate();
  if ((int) nx.size() == ndim) _nx = nx;
  if ((int) dx.size() == ndim) _dx = dx;
  if ((int) x0.size() == ndim) _x0 = x0;
}

Grid::Grid(const Grid &r)
  : AStringable(r)
{
  _recopy(r);
}

Grid& Grid::operator= (const Grid &r)
{
  AStringable::operator=(r);
  _recopy(r);
  return *this;
}

Grid::~Grid()
{
}

void Grid::resetFromSpaceDimension(int ndim)
{
  _nDim = ndim;
  _allocate();
}

int Grid::resetFromVector(const VectorInt& nx,
                          const VectorDouble& dx,
                          const VectorDouble& x0,
                          const VectorDouble& angles)
{
  _nDim = static_cast<int> (nx.size());
  _allocate();
  _nx = nx;
  for (int idim = 0; idim < _nDim; idim++)
  {
    if (nx[idim] < 0)
    {
      messerr("The number of grid mesh (%d) in direction (%d) may not be negative",
              nx[idim],idim+1);
      return 1;
    }
  }
  if (! x0.empty())
    _x0 = x0;
  else
    for (int idim = 0; idim < _nDim; idim++) _x0[idim] = 0.;
  if (! dx.empty())
  {
    _dx = dx;
    for (int idim = 0; idim < _nDim; idim++)
    {
      if (dx[idim] < 0.)
      {
        messerr("The mesh (%lf) in direction (%d) may not be negative",
            dx[idim],idim+1);
        return 1;
      }
    }
  }
  else
    for (int idim = 0; idim < _nDim; idim++) _dx[idim] = 1.;

  _rotation.setAngles(angles);
  return 0;
}

void Grid::resetFromGrid(Grid* grid)
{
  _nDim = grid->getNDim();
  _allocate();
  for (int idim=0; idim<_nDim; idim++)
  {
    _nx[idim] = grid->getNX(idim);
    _x0[idim] = grid->getX0(idim);
    _dx[idim] = grid->getDX(idim);
  }

  if (grid->isRotated())
  {
    setRotationByAngles(grid->getRotAngles());
  }
}

void Grid::setX0(int idim, double value)
{
  if (! _isSpaceDimensionValid(idim)) return;
  _x0[idim] = value;
}

void Grid::setDX(int idim, double value)
{
  if (! _isSpaceDimensionValid(idim)) return;
  if (value < 0) messageAbort("Argument 'dx' may not be negative");
  _dx[idim] = value;
}

void Grid::setNX(int idim, int value)
{
  if (! _isSpaceDimensionValid(idim)) return;
  if (value < 0) messageAbort("Argument 'nx' may not be negative");
  _nx[idim] = value;
}

void Grid::setRotationByMatrix(const MatrixSquareGeneral& rotmat)
{
  _rotation.resetFromSpaceDimension(_nDim);
  _rotation.setMatrixDirect(rotmat);
}

void Grid::setRotationByVector(const VectorDouble& rotmat)
{
  if (rotmat.empty()) return;
  _rotation.resetFromSpaceDimension(_nDim);
  _rotation.setMatrixDirectVec(rotmat);
}

void Grid::setRotationByAngles(const VectorDouble& angles)
{
  if (angles.empty()) return;
  _rotation.resetFromSpaceDimension(_nDim);
  _rotation.setAngles(angles);
}

/**
 * Define the rotation by the value of its first angle
 * @param angle Value of the first rotation angle
 */
void Grid::setRotationByAngle(double angle)
{
  _rotation.resetFromSpaceDimension(_nDim);
  VectorDouble angles(_nDim,0.);
  angles[0] = angle;
  _rotation.setAngles(angles);
}

double Grid::getX0(int idim) const
{
  if (! _isSpaceDimensionValid(idim)) return TEST;
  return _x0[idim];
}

double Grid::getDX(int idim) const
{
  if (! _isSpaceDimensionValid(idim)) return TEST;
  return _dx[idim];
}

int Grid::getNX(int idim) const
{
  if (! _isSpaceDimensionValid(idim)) return ITEST;
  return _nx[idim];
}

int Grid::getNTotal() const
{
  if (_nDim <= 0) return 0;
  int ntotal = 1;
  for (int idim=0; idim<_nDim; idim++) 
    ntotal *= _nx[idim];
  return ntotal;
}

double Grid::getVolume(bool flagCell) const
{
  double volume = 1;
  for (int idim = 0; idim < _nDim; idim++)
    volume *= getExtend(idim, flagCell);
  return volume;
}

double Grid::getExtend(int idim, bool flagCell) const
{
  double ext;
  if (flagCell)
    ext = _nx[idim] * _dx[idim];
  else
    ext = (_nx[idim] - 1) * _dx[idim];
  return ext;
}

VectorDouble Grid::getExtends(bool flagCell) const
{
  VectorDouble ext(_nDim);
  for (int idim = 0; idim < _nDim; idim++)
      ext[idim] = getExtend(idim, flagCell);
  return ext;
}

double Grid::getCellSize() const
{
  double size = 1.;
  for (int idim=0; idim<_nDim; idim++)
    size *= _dx[idim];
  return size;
}

String Grid::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  if (_nDim <= 0) return sstr.str();

  sstr << toTitle(1,"Grid characteristics:");
  sstr << "Origin : ";
  for (int idim=0; idim<_nDim; idim++)
    sstr << toDouble(_x0[idim]);
  sstr << std::endl;

  sstr << "Mesh   : ";
  for (int idim=0; idim<_nDim; idim++)
    sstr << toDouble(_dx[idim]);
  sstr << std::endl;

  sstr << "Number : ";
  for (int idim=0; idim<_nDim; idim++)
    sstr << toInt(_nx[idim]);
  sstr << std::endl;

  sstr << _rotation.toString(strfmt);

  return sstr.str();
}

double Grid::getCoordinate(int rank, int idim0, bool flag_rotate) const
{
  /* Convert a sample number into grid indices */

  rankToIndice(rank, _iwork0);

  /* Calculate the coordinates in the grid system */

  for (int idim = 0; idim < _nDim; idim++)
    _work1[idim] = _iwork0[idim] * _dx[idim];

  /* Process the grid rotation (if any) */

  if (flag_rotate)
  {
    _rotation.rotateDirect(_work1, _work2);
    return (_work2[idim0] + _x0[idim0]);
  }
  return (_work1[idim0] + _x0[idim0]);
}

/**
 * Returns the coordinates of a grid node, defined by its indices
 * @param indice       Vector of indices defining the target grid node
 * @param flag_rotate  True if the grid rotation must be taken into account
 * @param shift        Vector of shifts (dimension: ndim)
 *                     0 : no shift; -1 : minus half a cell-width; +1 plus half a cell-width
 * @param dxsPerCell   Vector of variable grid meshes (optional)
 * @return
 */
VectorDouble Grid::getCoordinatesByIndice(const VectorInt &indice,
                                          bool flag_rotate,
                                          const VectorInt& shift,
                                          const VectorDouble& dxsPerCell) const
{
  /* Calculate the coordinates in the grid system */

  for (int idim = 0; idim < _nDim; idim++)
  {
    _work1[idim] = indice[idim] * _dx[idim];
    if (! shift.empty())
    {
      double ext = (dxsPerCell.empty()) ? _dx[idim] : dxsPerCell[idim];
      _work1[idim] += shift[idim] * ext / 2.;
    }
  }

  /* Process the grid rotation (if any) */

  if (flag_rotate)
  {
    _rotation.rotateDirect(_work1, _work2);
    for (int idim = 0; idim < _nDim; idim++) _work2[idim] += _x0[idim];
    return _work2;
  }
  for (int idim = 0; idim < _nDim; idim++) _work1[idim] += _x0[idim];
  return _work1;
}

/**
 * Returns the coordinates of a Grid corner
 * @param icorner Vector specifying the corner (0: minimum; 1: maximum). (Dimension: ndim)
 * @return The coordinates of a corner
 */
VectorDouble Grid::getCoordinatesByCorner(const VectorInt& icorner) const
{
  VH::fill(_iwork0, 0);
  for (int idim = 0; idim < _nDim; idim++)
    if (icorner[idim] > 0) _iwork0[idim] = _nx[idim]-1;
  return getCoordinatesByIndice(_iwork0);
}

/**
 * Returns the coordinates of a Grid cell corner
 * @param node Rank of the Target cell
 * @param shift        Vector of shifts (dimension: ndim)
 *                     0 : no shift; -1 : minus half a cell-width; +1 plus half a cell-width
 * @param dxsPerCell Vector of variable mesh extensions at target cell
 * @return The coordinates of a cell corner (possibly shifted)
 */
VectorDouble Grid::getCellCoordinatesByCorner(int node,
                                              const VectorInt& shift,
                                              const VectorDouble& dxsPerCell) const
{
  rankToIndice(node, _iwork0);
  return getCoordinatesByIndice(_iwork0, true, shift, dxsPerCell);
}

/**
 * Return the Vector of coordinates for a given grid node
 * @param rank Rank of the target grid node
 * @param flag_rotate TRUE: perform the rotation; FALSE: skip rotation
 * @return Vector of coordinates
 */
VectorDouble Grid::getCoordinatesByRank(int rank, bool flag_rotate) const
{
  /* Convert a sample number into grid indices */

  rankToIndice(rank, _iwork0);

  /* Calculate the coordinates in the grid system */

  for (int idim = 0; idim < _nDim; idim++)
    _work1[idim] = _iwork0[idim] * _dx[idim];

  /* Process the grid rotation (if any) */

  if (flag_rotate)
  {
    _rotation.rotateDirect(_work1, _work2);
    for (int idim = 0; idim < _nDim; idim++) _work2[idim] += _x0[idim];
    return _work2;
  }

  /* Shift for the origin */

  for (int idim = 0; idim < _nDim; idim++)
    _work1[idim] += _x0[idim];
  return _work1;
}

double Grid::indiceToCoordinate(int idim0,
                                const constvectint indice,
                                const constvect percent,
                                bool flag_rotate) const
{
  /* Calculate the coordinates in the grid system */

  for (int idim=0; idim<_nDim; idim++)
  {
    _work1[idim] = (double) indice[idim];
    if (! percent.empty()) _work1[idim] += percent[idim];
    _work1[idim] *= _dx[idim];
  }

  if (flag_rotate)
  {
    _rotation.rotateDirect(_work1, _work2);
    return (_work2[idim0] + _x0[idim0]);
  }
  return (_work1[idim0] + _x0[idim0]);
}

VectorDouble Grid::indicesToCoordinate(const VectorInt& indice,
                                       const VectorDouble& percent) const
{
  indicesToCoordinateInPlace(indice, _work1, percent);
  return _work1;
}

void Grid::indicesToCoordinateInPlace(const VectorInt& indice,
                                      VectorDouble& coor,
                                      const VectorDouble& percent,
                                      bool flag_rotate) const
{
  if ((int)coor.size() < _nDim)
  {
    messerr("Argument coor should have the correct size. Output argument 'coor' not modified.");
    return;
  }

  /* Calculate the coordinates in the grid system */

  for (int idim=0; idim<_nDim; idim++)
  {
    _work1[idim] = indice[idim];
    if (! percent.empty()) _work1[idim] += percent[idim];
    _work1[idim] *= _dx[idim];
  }

  if (flag_rotate)
  {
    _rotation.rotateDirect(_work1,_work2);
    for (int idim = 0; idim < _nDim; idim++)
      coor[idim] = _work2[idim] + _x0[idim];
  }
  else
  {
    for (int idim = 0; idim < _nDim; idim++)
      coor[idim] = _work1[idim] + _x0[idim];
  }
}

double Grid::rankToCoordinate(int idim0, int rank, const VectorDouble& percent) const
{
  rankToIndice(rank, _iwork0);
  return indiceToCoordinate(idim0, _iwork0, percent);
}

VectorDouble Grid::rankToCoordinates(int rank, const VectorDouble& percent) const
{
  rankToIndice(rank, _iwork0);
  return indicesToCoordinate(_iwork0,percent);
}

void Grid::rankToCoordinatesInPlace(int rank, VectorDouble& coor, const VectorDouble& percent) const
{
  rankToIndice(rank, _iwork0);
  return indicesToCoordinateInPlace(_iwork0, coor, percent);
}

int Grid::indiceToRank(const constvectint indice) const
{
  int ival = indice[_nDim-1];
  if (ival < 0 || ival >= _nx[_nDim-1]) return(-1);
  for (int idim=_nDim-2; idim>=0; idim--)
  {
    if (indice[idim] < 0 || indice[idim] >= _nx[idim]) return(-1);
    ival = ival * _nx[idim] + indice[idim];
  }
  return ival;
}

/**
 *
 * @param rank Rank of the Node (in the meshing)
 * @param indices Indices of the node in the grid system
 * @param minusOne Consider that the number of cells in each direction
 * should be reduced by one.
 *
 * \remarks: The number of nodes in the grid per direction
 * \remarks: must be adapted (subtracting 1) due to interval.
 */
void Grid::rankToIndice(int rank, vectint indices, bool minusOne) const
{
  int minus = (minusOne) ? 1 : 0;

  const int* nxadd = _nx.data(); // for optimization, use address rather than []
  int nval         = 1;
  for (int idim = 0; idim < _nDim; idim++)
    nval *= (*(nxadd + idim) - minus);

  int newind;
  for (int idim = _nDim - 1; idim >= 0; idim--)
  {
    nval /= (*(nxadd + idim) - minus);
    newind           = rank / nval;
    indices[idim] = newind;
    rank -= newind * nval;
  }
}

VectorInt Grid::coordinateToIndices(const VectorDouble &coor,
                                    bool centered,
                                    double eps) const
{
  if (coordinateToIndicesInPlace(coor, _iwork0, centered, eps)) return VectorInt();
  return _iwork0;
}

/**
 * Find the grid node to which the current sample is assigned
 * @param coor     Sample coordinates
 * @param indice   Indices of the assigned grid node
 * @param centered True for grid cell centered
 * @param eps      Epsilon to over-pass roundoff problem
 * @return Error return code
 */
int Grid::coordinateToIndicesInPlace(const VectorDouble &coor,
                                     VectorInt &indice,
                                     bool centered,
                                     double eps) const
{
  if ((int)indice.size() != _nDim)
  {
    messerr("Argument 'indice' should have the correct size. Output argument 'indice' not modified.");
    return -1;
  }

  // Check if all coordinates are defined

  for (int idim = 0; idim < _nDim; idim++)
    if (FFFF(coor[idim])) return -1;

  // Shift by the origin

  for (int idim = 0; idim < _nDim; idim++)
    _work1[idim] = coor[idim] - _x0[idim];

  // Perform the Inverse rotation

  _rotation.rotateInverse(_work1, _work2);

  // Calculate the indices

  bool outside = false;
  for (int idim = 0; idim < _nDim; idim++)
  {
    int ix;
    if (centered)
      ix = (int) floor(_work2[idim] / _dx[idim] + 0.5 + eps);
    else
      ix = (int) floor(_work2[idim] / _dx[idim] + eps);
    indice[idim] = ix;
    if (ix < 0 || ix >= _nx[idim]) outside = true;
  }
  return (int) outside;
}

int Grid::coordinateToRank(const VectorDouble& coor, bool centered, double eps) const
{
  if (coordinateToIndicesInPlace(coor,_iwork0,centered,eps)) return -1;
  return indiceToRank(_iwork0);
}

VectorInt Grid::getCenterIndices() const
{
  for (int idim = 0; idim < _nDim; idim++)
    _iwork0[idim] = (_nx[idim] - 1) / 2;
  return _iwork0;
}

bool Grid::_isSpaceDimensionValid(int idim) const
{
  return checkArg("Argument 'idim' is invalid", idim, _nDim);
}

void Grid::_allocate(void)
{
  _nx.resize(_nDim);
  for (int i=0; i<_nDim; i++) _nx[i] = 1;
  _x0.resize(_nDim);
  for (int i=0; i<_nDim; i++) _x0[i] = 0.;
  _dx.resize(_nDim);
  for (int i=0; i<_nDim; i++) _dx[i] = 0.;

  _rotation.resetFromSpaceDimension(_nDim);

  _iwork0.resize(_nDim);
  _work1.resize(_nDim);
  _work2.resize(_nDim);
}

void Grid::_recopy(const Grid &r)
{
  _nDim = r._nDim;
  _allocate();

  _nx = r._nx;
  _x0 = r._x0;
  _dx = r._dx;

  // Copy the rotation

  _rotation = r._rotation;

  _iter = r._iter;
  _nprod = r._nprod;
  _counts = r._counts;
  _order = r._order;
  _indices = r._indices;

  // For working array, simply dimension but do not loose time in copying the contents

  _iwork0 = VectorInt(_nDim);
  _work1 = VectorDouble(_nDim);
  _work2 = VectorDouble(_nDim);
}

/**
 * Copy some parameters from Gridaux
 * @param mode     Type of parameters to be copied
 *                 1 : Array of Grid Number of meshes
 *                 2 : Array of Grid origin
 *                 3 : Array of Grid Meshes
 *                 4 : Rotation
 * @param gridaux  Source Grid structure
 */
void Grid::copyParams(int mode, const Grid& gridaux)
{
  if (gridaux.getNDim() != _nDim)  return;

  /* Dispatch */

  switch (mode)
  {
    case 1: /* Copy NX */
      _nx = gridaux.getNXs();
      break;

    case 2: /* Copy X0 */
      _x0 = gridaux.getX0s();
      break;

    case 3: /* Copy DX */
      _dx = gridaux.getDXs();
      break;

    case 4: /* Rotation matrices */
      setRotationByAngles(gridaux.getRotAngles());
  }
}

/**
 * Check that the current grid match the one provided as argument
 * up to their common Space Dimension
 * @param grid Target grid to be checked against the current one
 * @return True if the grid match
 */
bool Grid::isSame(const Grid& grid) const
{
  int ndim = MIN(_nDim, grid.getNDim());

  /* Compare the grid parameters */

  for (int idim = 0; idim < ndim; idim++)
  {
    if (_nx[idim] != grid.getNX(idim)) return 0;
    if (_dx[idim] != grid.getDX(idim)) return 0;
    if (_x0[idim] != grid.getX0(idim)) return 0;
  }

  /* Compare the rotations */

  if (isRotated() != grid.isRotated()) return 0;
  if (isRotated())
  {
    for (int idim = 0; idim < ndim; idim++)
      if (getRotAngle(idim) != grid.getRotAngle(idim)) return 0;
  }
  return 1;
}

bool Grid::isSameMesh(const Grid& grid) const
{
  int ndim = MIN(_nDim, grid.getNDim());

  /* Compare the grid parameters */

  for (int idim = 0; idim < ndim; idim++)
  {
    if (_dx[idim] != grid.getDX(idim)) return 0;
  }
  return 1;
}

/**
 * Returns a vector with the coordinates along one axis. This is needed
 * for the label of Grid representation
 * Warning: Not considering any possible rotation.
 * @param idim Index of the Space Dimension
 * @return
 */
VectorDouble Grid::getAxis(int idim) const
{
  VectorDouble vect;
  if (idim < 0 || idim >= getNDim()) return (vect);

  int nvect     = getNX(idim);
  double origin = getX0(idim);
  double pas    = getDX(idim);
  vect.resize(nvect);

  for (int i = 0; i < nvect; i++)
    vect[i] = origin + i * pas;
  return vect;
}

/**
 * Initialize an iterator on the grid
 * @param order Array giving the order of the Space Dimensions when iterating
 */
void Grid::iteratorInit(const VectorInt& order)
{
  // Initiate the Iterator
  _iter = 0;

  // Define the number of positions per space index
  _counts = _nx;

  // Define the order
  if (order.empty() || _nDim != (int) order.size())
  {
    _order.resize(_nDim,0);
    for (int idim = 0; idim < _nDim; idim++)
      _order[idim] = idim;
  }
  else
  {
    // Check that the order array provided is valid
    for (int idim = 0; idim < _nDim; idim++)
    {
      bool found = false;
      for (int jdim = 0; jdim < _nDim; jdim++)
      {
        int rank = ABS(order[jdim]) - 1;
        if (rank == idim) found = true;
      }
      if (! found)
      {
        messerr("When provided, 'order' should contain all Space dimensions. Iterator cancelled.");
        _iter = 0;
        _counts = VectorInt();
        _order = VectorInt();
        return;
      }
    }
    _order = order;
  }

  // Count the total number of iterations
  _nprod = 1;
  for (int idim = 0; idim < _nDim; idim++)
    _nprod *= _counts[idim];
}

/**
 * Return the vector of grid indices for each iteration
 * @return
 */
VectorInt Grid::iteratorNext(void)
{
  int idim;
  int iech = _iter;
  int nval = _nprod;

  VectorInt indices(_nDim);
  for (int jdim = _nDim - 1; jdim >= 0; jdim--)
  {
    int order = _order[jdim];
    idim = ABS(order);
    nval /= _counts[idim];
    int divid = iech / nval;
    indices[idim] = divid;
    iech -= divid * nval;
  }

  // Increment the iterator
  if (_iter < _nprod - 1) _iter++;
  return indices;
}

bool Grid::empty() const
{
  bool empty = _nDim <= 0 || _nx.empty();
  return empty;
}

/****************************************************************************/
/*!
 **  Returns the characteristics of a dilated grid
 **
 ** \param[in]  mode   1 for extending; -1 for compressing
 ** \param[in]  nshift Array of shifts
 **
 ** \param[out] nx    Array of number of grid meshes
 ** \param[out] dx    Array of grid meshes
 ** \param[out] x0    Array of grid origins
 **
 *****************************************************************************/
void Grid::dilate(int mode,
                  const VectorInt& nshift,
                  VectorInt& nx,
                  VectorDouble& dx,
                  VectorDouble& x0) const
{
  if (mode != 1 && mode != -1) return;

  /* Get the number of grid nodes */

  for (int idim = 0; idim < _nDim; idim++)
  {
    nx[idim] = getNX(idim) + 2 * mode * nshift[idim];
    if (nx[idim] <= 0) return;
    dx[idim] = getDX(idim);
  }

  /* Get the lower left corner of the small grid */

  for (int idim = 0; idim < _nDim; idim++)
    _iwork0[idim] = -mode * nshift[idim];
  indicesToCoordinate(_iwork0, _work1);

  /* Calculate the center of the lower left cell */

  for (int idim = 0; idim < _nDim; idim++)
    x0[idim] = _work1[idim];
}

/****************************************************************************/
/*!
 **  Returns the characteristics of a multiple grid
 **
 ** \param[in]  nmult Array of multiplicity coefficients
 ** \param[in]  flagCell true for cell matching; 0 for point matching
 **
 ** \param[out] nx    Array of number of grid meshes
 ** \param[out] dx    Array of grid meshes
 ** \param[out] x0    Array of grid origins
 **
 *****************************************************************************/
void Grid::multiple(const VectorInt &nmult,
                    bool flagCell,
                    VectorInt &nx,
                    VectorDouble &dx,
                    VectorDouble &x0) const
{
  VectorDouble perc(_nDim);
  VectorDouble coor1(_nDim); // cannot use _work1 and _work2 (as already used inside)
  VectorDouble coor2(_nDim);

  /* Get the number of grid nodes */

  for (int idim = 0; idim < _nDim; idim++)
  {
    double value = (double) getNX(idim);
    if (flagCell)
      nx[idim] = (int) floor(value / (double) nmult[idim]);
    else
      nx[idim] = 1 + (int) floor((value - 1.) / (double) nmult[idim]);
  }

  /* Get the new grid meshes */

  for (int idim = 0; idim < _nDim; idim++)
    dx[idim] = getDX(idim) * nmult[idim];

  /* Get the lower left corner of the small grid */

  for (int idim = 0; idim < _nDim; idim++) _iwork0[idim] = 0;
  for (int idim = 0; idim < _nDim; idim++) perc[idim] = -0.5;
  indicesToCoordinateInPlace(_iwork0, coor1, perc);
  for (int idim = 0; idim < _nDim; idim++) perc[idim] = 0.5;
  indicesToCoordinateInPlace(_iwork0, coor2, perc);

  /* Calculate the center of the lower left cell */

  for (int idim = 0; idim < _nDim; idim++)
  {
    double delta = (coor2[idim] - coor1[idim]) / 2.;
    if (flagCell)
      x0[idim] = coor1[idim] + delta * (double) nmult[idim];
    else
      x0[idim] = getX0(idim);
  }
}

/****************************************************************************/
/*!
 **  Returns the characteristics of a divider grid
 **
 ** \param[in]  nmult Array of subdivision coefficients
 ** \param[in]  flagCell true for cell matching; 0 for point matching
 **
 ** \param[out] nx    Array of number of grid meshes
 ** \param[out] dx    Array of grid meshes
 ** \param[out] x0    Array of grid origins
 **
 *****************************************************************************/
void Grid::divider(const VectorInt &nmult,
                   bool flagCell,
                   VectorInt &nx,
                   VectorDouble &dx,
                   VectorDouble &x0) const
{
  VectorDouble perc(_nDim);
  VectorDouble coor1(_nDim); // cannot use _work1 and _work2 (as already used inside)
  VectorDouble coor2(_nDim);

  /* Get the number of grid nodes */

  for (int idim = 0; idim < _nDim; idim++)
  {
    if (flagCell)
      nx[idim] = getNX(idim) * nmult[idim];
    else
      nx[idim] = 1 + (getNX(idim) - 1) * nmult[idim];
  }

  /* Get the new grid meshes */

  for (int idim = 0; idim < _nDim; idim++)
    dx[idim] = getDX(idim) / ((double) nmult[idim]);

  /* Get the lower left corner of the small grid */

  for (int idim = 0; idim < _nDim; idim++) _iwork0[idim] = 0;
  for (int idim = 0; idim < _nDim; idim++) perc[idim] = -0.5;
  indicesToCoordinateInPlace(_iwork0, coor1, perc);
  for (int idim = 0; idim < _nDim; idim++) perc[idim] = 0.5;
  indicesToCoordinateInPlace(_iwork0, coor2, perc);

  /* Calculate the center of the lower left cell */

  for (int idim = 0; idim < _nDim; idim++)
  {
    double delta = (coor2[idim] - coor1[idim]) / 2.;
    if (flagCell)
      x0[idim] = coor1[idim] + delta / (double) nmult[idim];
    else
      x0[idim] = getX0(idim);
  }
}

/****************************************************************************/
/*!
 **  Return the index of a sample when calculated from mirroring within
 **  an array whose indices vary between 0 and nx-1
 **
 ** \return Rank of the restrained cell
 **
 ** \param[in]  idim      Rank of the space dimension
 ** \param[in]  ix        Rank of the cell to be restrained
 **
 *****************************************************************************/
int Grid::getMirrorIndex(int idim, int ix) const
{
  return generateMirrorIndex(_nx[idim], ix);
}

/****************************************************************************/
/*!
 **  Returns an array giving the ranks of the nodes (according to user's order)
 **  coded with standard order (according to gstlearn internal order)
 **
 ** \return Array of indices
 **
 ** \param[in]  nx      Array giving the number of cells per direction
 ** \param[in]  string  String describing the sorting order
 ** \param[in]  startFromZero True if numbering must start from 0 (1 otherwise)
 ** \param[in]  invert  Way to use the resulting array (see remark)
 ** \param[in]  verbose Verbose flag
 **
 ** \remark Example of string: "+x2-x1"
 **
 ** \remark if 'rank' designates the resulting vector of indices
 ** \remark invert=True:
 ** \remark   rank[i] is the location of element 'i' of the user's array
 ** \remark   within a regular grid of gstlearn
 ** \remark invert=False:
 ** \remark   rank[i] is the rank of the element of the user's array
 ** \remark   in position 'i' of the regular grid of gstlearn
 **
 *****************************************************************************/
VectorInt Grid::gridIndices(const VectorInt &nx,
                            const String &string,
                            bool startFromZero,
                            bool invert,
                            bool verbose)
{
  int ndim = (int) nx.size();
  int ncell = VH::product(nx);

  // Decode the string

  VectorInt order = decodeGridSorting(string, nx, verbose);
  if (order.empty()) return VectorInt();

  // Initialization of the recursion structure

  DimLoop dlp;
  dlp.curech = 0;
  dlp.ndim = ndim;
  dlp.nx = nx;
  dlp.order = order;
  dlp.indg = VectorInt(ndim,0);
  dlp.tab = VectorInt(ncell);

  // Recursion

  _dimensionRecursion(ndim-1, verbose, dlp);

  // Invert order

  VectorInt tab2(dlp.tab);
  if (invert)
  {
    VectorInt ind(VH::sequence(ncell));
    VH::arrangeInPlace(0, ind, dlp.tab, true, ncell);
    for (int i = 0; i < ncell; i++)
    {
      tab2[i] = dlp.tab[ind[i]];
    }
  }

  // Change the starting index
  if (startFromZero)
  {
    for (int i = 0; i < ncell; i++)
    {
      tab2[i] -= 1;
    }
  }
  return tab2;
}

VectorInt Grid::generateGridIndices(const String &string,
                                    bool startFromZero,
                                    bool invert,
                                    bool verbose) const
{
  return gridIndices(getNXs(), string, startFromZero, invert, verbose);
}

/****************************************************************************/
/*!
 **  Return the index of a sample when calculated from mirroring within
 **  an array whose indices vary between 0 and nx-1
 **
 ** \return Rank of the restrained cell
 **
 ** \param[in]  nx        Number of cells
 ** \param[in]  ix        Rank of the cell to be restrained
 **
 *****************************************************************************/
int Grid::generateMirrorIndex(int nx, int ix)
{
  int nmax = nx - 1;
  while (ix < 0 || ix >= nx)
  {
    if (ix < 0)
    {
      ix = -ix;
    }
    else if (ix > nmax)
    {
      ix = 2 * nmax - ix;
    }
  }
  return (ix);
}

/**
 * Check if a sample belongs to a Grid Cell
 * @param coor       Sample coordinates
 * @param center     Coordinates of the grid node center
 * @param dxsPerCell When defined, vector of cell extension; otherwise use dx
 * @return Error return code
 *
 * @remark Samples located exactly on the edge are considered as INSIDE
 */
bool Grid::sampleBelongsToCell(constvect coor,
                               constvect center,
                               const VectorDouble& dxsPerCell) const
{
  if (_rotation.isRotated())
  {
    // Convert Grid Node center into Grid coordinates
    VectorDouble _work3 = _work1;
    for (int idim = 0; idim < _nDim; idim++)
      _work1[idim] = center[idim] - _x0[idim];
    _rotation.rotateInverse(_work1, _work3);

    // Convert the coordinates of sample in Grid coordinates
    for (int idim = 0; idim < _nDim; idim++)
      _work1[idim] = coor[idim] - _x0[idim];
    _rotation.rotateInverse(_work1, _work2);

    // Calculate the departure between sample and grid center
    for (int idim = 0; idim < _nDim; idim++)
    {
      double dxloc = (dxsPerCell.empty()) ? _dx[idim] : dxsPerCell[idim];
      double delta = _work2[idim] - _work3[idim];
      if (ABS(delta) > dxloc / 2.) return false;
    }
  }
  else
  {
    // Calculate the departure between sample and grid center
    for (int idim = 0; idim < _nDim; idim++)
    {
      double dxloc = (dxsPerCell.empty()) ? _dx[idim] : dxsPerCell[idim];
      double delta = center[idim] - coor[idim];
      if (ABS(delta) > dxloc / 2.) return false;
    }
  }
  return true;
}

bool Grid::sampleBelongsToCell(const VectorDouble& coor,
                               const VectorDouble& center,
                               const VectorDouble& dxsPerCell) const
{
 return sampleBelongsToCell(constvect(coor),constvect(center),dxsPerCell);
}

/**
 * Check if a sample belongs to a Grid Cell
 * @param coor       Sample coordinates (can be lower space dimension than the current Grid)
 * @param rank       Rank of the Grid cell
 * @param dxsPerCell When defined, vector of cell extension; otherwise use dx
 * @return Error return code
 *
 * @remark Samples located exactly on the edge are considered as INSIDE
 */
bool Grid::sampleBelongsToCell(const VectorDouble &coor,
                               int rank,
                               const VectorDouble &dxsPerCell) const
{
  // Identify the coordinates of the center of the grid cell, referred by its 'rank' and
  // convert into Grid coordinates
  VectorDouble center = rankToCoordinates(rank);

  // Complement 'coor' to the grid space dimension
  int ndim_coor = (int) coor.size();
  VectorDouble coor_loc;
  if (ndim_coor == _nDim)
    coor_loc = coor;
  else
  {
    coor_loc = center;
    for (int idim = 0; idim < ndim_coor; idim++)
      coor_loc[idim] = coor[idim];
  }

  if (_rotation.isRotated())
  {
    VectorDouble _work3 = _work1;
    for (int idim = 0; idim < _nDim; idim++)
      _work1[idim] = center[idim] - _x0[idim];
    _rotation.rotateInverse(_work1, _work3);

    // Convert the coordinates of sample in Grid coordinates
    for (int idim = 0; idim < _nDim; idim++)
      _work1[idim] = coor_loc[idim] - _x0[idim];
    _rotation.rotateInverse(_work1, _work2);

    // Calculate the departure between sample and grid center
    for (int idim = 0; idim < MIN(ndim_coor, _nDim); idim++)
    {
      double dxloc = (dxsPerCell.empty()) ? _dx[idim] : dxsPerCell[idim];
      double delta = _work2[idim] - _work3[idim];
      if (ABS(delta) > dxloc / 2.) return false;
    }
  }
  else
  {
    // Calculate the departure between sample and grid center
     for (int idim = 0; idim < MIN(ndim_coor, _nDim); idim++)
     {
       double dxloc = (dxsPerCell.empty()) ? _dx[idim] : dxsPerCell[idim];
       double delta = center[idim] - coor_loc[idim];
       if (ABS(delta) > dxloc / 2.) return false;
     }
  }
  return true;
}
