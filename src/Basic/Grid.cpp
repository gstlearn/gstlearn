/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Basic/Rotation.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Grid.hpp"
#include "Basic/Grid.hpp"

#include <math.h>

Grid::Grid(int ndim,
             const VectorInt& nx,
             const VectorDouble& x0,
             const VectorDouble& dx)
  : AStringable(),
    _nDim(ndim)
  , _nx(nx)
  , _x0(x0)
  , _dx(dx)
  , _rotation()
  , _iter(0)
  , _nprod(0)
  , _counts()
  , _order()
  , _indices()
{
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
      messerr("The number of grid mesh (%d) is direction (%d) may not be negative",
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

void Grid::setX0(int idim,
                  double value)
{
  if (! _isSpaceDimensionValid(idim)) return;
  _x0[idim] = value;
}

void Grid::setDX(int idim,
                  double value)
{
  if (! _isSpaceDimensionValid(idim)) return;
  if (value < 0) messageAbort("Argument 'dx' may not be negative");
  _dx[idim] = value;
}

void Grid::setNX(int idim,
                  int value)
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
  _rotation.setMatrixDirectByVector(rotmat);
}

void Grid::setRotationByAngles(const VectorDouble angles)
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
  int ntotal = 1;
  for (int idim=0; idim<_nDim; idim++) 
    ntotal *= _nx[idim];
  return ntotal;
}

double Grid::getVolume(bool flag_cell) const
{
  double volume = 1;
  for (int idim = 0; idim < _nDim; idim++)
    volume *= getExtend(idim, flag_cell);
  return volume;
}

double Grid::getExtend(int idim, bool flag_cell) const
{
  double ext;
  if (flag_cell)
    ext = _nx[idim] * _dx[idim];
  else
    ext = (_nx[idim] - 1) * _dx[idim];
  return ext;
}

VectorDouble Grid::getExtends(bool flag_cell) const
{
  VectorDouble ext(_nDim);
  for (int idim = 0; idim < _nDim; idim++)
      ext[idim] = getExtend(idim, flag_cell);
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
  int ndim = _nDim;

  sstr << toTitle(1,"Grid characteristics:");
  sstr << "Origin : ";
  for (int idim=0; idim<ndim; idim++)
    sstr << toDouble(_x0[idim]);
  sstr << std::endl;

  sstr << "Mesh   : ";
  for (int idim=0; idim<ndim; idim++)
    sstr << toDouble(_dx[idim]);
  sstr << std::endl;

  sstr << "Number : ";
  for (int idim=0; idim<ndim; idim++)
    sstr << toInt(_nx[idim]);
  sstr << std::endl;

  sstr << _rotation.toString(strfmt);

  return sstr.str();
}

double Grid::getCoordinate(int rank, int idim0, bool flag_rotate) const
{
  int ndim = getNDim();
  VectorInt   iwork(ndim);
  VectorDouble work1(ndim);
  VectorDouble work2(ndim);

  /* Convert a sample number into grid indices */

  rankToIndice(rank, iwork);

  /* Calculate the coordinates in the grid system */

  for (int idim = 0; idim < ndim; idim++)
    work1[idim] = iwork[idim] * _dx[idim];

  /* Process the grid rotation (if any) */

  if (flag_rotate)
    _rotation.rotateInverse(work1,work2);
  else
    work2 = work1;

  return (work2[idim0] + _x0[idim0]);
}

VectorDouble Grid::getCoordinatesByIndice(const VectorInt& indice, bool flag_rotate) const
{
  int ndim = getNDim();
  VectorDouble work1(ndim);
  VectorDouble work2(ndim);

  /* Calculate the coordinates in the grid system */

  for (int idim = 0; idim < ndim; idim++)
    work1[idim] = indice[idim] * _dx[idim];

  /* Process the grid rotation (if any) */

  if (flag_rotate)
    _rotation.rotateInverse(work1,work2);
  else
    work2 = work1;

  for (int idim = 0; idim < ndim; idim++)
    work2[idim] += _x0[idim];

  return work2;
}

/**
 * Returns the coordinates of a Grid corner
 * @param icorner Vector specifying the corner (0: minimum; 1: maximum). (Dimension: ndim)
 * @return The coordinates of a corner
 */
VectorDouble Grid::getCoordinatesByCorner(const VectorInt& icorner) const
{
  VectorInt indice(_nDim,0);
  for (int idim = 0; idim < _nDim; idim++)
    if (icorner[idim] > 0) indice[idim] = _nx[idim]-1;
  return getCoordinatesByIndice(indice);
}

/**
 * Return the Vector of coordinates for a given grid node
 * @param rank Rank of the target grid node
 * @param flag_rotate TRUE: perform the roataion; FALSE: skip rotation
 * @return Vector of coordinates
 */
VectorDouble Grid::getCoordinatesByRank(int rank, bool flag_rotate) const
{
  int ndim = getNDim();
  VectorInt    iwork(ndim);
  VectorDouble work1(ndim);
  VectorDouble work2(ndim);

  /* Convert a sample number into grid indices */

  rankToIndice(rank, iwork);

  /* Calculate the coordinates in the grid system */

  for (int idim = 0; idim < ndim; idim++)
    work1[idim] = iwork[idim] * _dx[idim];

  /* Process the grid rotation (if any) */

  if (flag_rotate)
    _rotation.rotateInverse(work1,work2);
  else
    work2 = work1;

  for (int idim = 0; idim < ndim; idim++)
    work2[idim] += _x0[idim];

  return work2;
}

double Grid::indiceToCoordinate(int idim0,
                                 const VectorInt& indice,
                                 const VectorDouble& percent) const
{
  VectorDouble work1(_nDim);
  VectorDouble work2(_nDim);

  /* Calculate the coordinates in the grid system */

  for (int idim=0; idim<_nDim; idim++)
  {
    work1[idim] = (double) indice[idim];
    if (! percent.empty()) work1[idim] += percent[idim];
    work1[idim] *= _dx[idim];
  }

  /* Process the grid rotation (if any) */

  _rotation.rotateInverse(work1,work2);

  // Shift by the origin of the grid

  return (work2[idim0] + _x0[idim0]);
}

VectorDouble Grid::indicesToCoordinate(const VectorInt& indice,
                                       const VectorDouble& percent) const
{
  VectorDouble vect(_nDim);
  indicesToCoordinateInPlace(indice, vect, percent);
  return vect;
}

void Grid::indicesToCoordinateInPlace(const VectorInt& indice,
                                      VectorDouble& coor,
                                      const VectorDouble& percent) const
{
  VectorDouble work1(_nDim);
  VectorDouble work2(_nDim);

  /* Calculate the coordinates in the grid system */

  for (int idim=0; idim<_nDim; idim++)
  {
    work1[idim] = indice[idim];
    if (! percent.empty()) work1[idim] += percent[idim];
    work1[idim] *= _dx[idim];
  }

  /* Process the grid rotation (if any) */

  _rotation.rotateInverse(work1,work2);

  // Returning vector

  for (int idim = 0; idim < _nDim; idim++)
    coor[idim] = work2[idim] + _x0[idim];
}

double Grid::rankToCoordinate(int idim0, int rank, const VectorDouble& percent) const
{
  VectorInt indice;
  rankToIndice(rank, indice);
  return indiceToCoordinate(idim0, indice, percent);
}

VectorDouble Grid::rankToCoordinates(int rank, const VectorDouble& percent) const
{
  VectorInt indice;
  rankToIndice(rank, indice);
  return indicesToCoordinate(indice,percent);
}

void Grid::rankToCoordinatesInPlace(int rank, VectorDouble& coor, const VectorDouble& percent) const
{
  VectorInt indice(_nDim);
  rankToIndice(rank, indice);
  return indicesToCoordinateInPlace(indice, coor, percent);
}

int Grid::indiceToRank(const VectorInt& indice) const
{
  int ndim = _nDim;
  int ival = indice[ndim-1];
  if (ival < 0 || ival >= _nx[ndim-1]) return(-1);
  for (int idim=ndim-2; idim>=0; idim--)
  {
    if (indice[idim] < 0 || indice[idim] >= _nx[idim]) return(-1);
    ival = ival * _nx[idim] + indice[idim];
  }
  return ival;
}

/**
 *
 * @param rank Rank of the Node (in the meshing)
 * @param indice Indices of the node in the grid system
 * @param minusOne Consider that the number of cells in each direction
 * should be reduced by one.
 *
 * \remarks: The number of nodes in the grid per direction
 * \remarks: must be adapted (subtracting 1) due to interval.
 */
void Grid::rankToIndice(int  rank,
                         VectorInt& indice,
                         bool minusOne) const
{
  int ndim = _nDim;
  int minus = (minusOne) ? 1 : 0;
  int nval = 1;
  for (int idim=0; idim<ndim; idim++) nval *= (_nx[idim] - minus);

  for (int idim=ndim-1; idim>=0; idim--)
  {
    nval /= (_nx[idim] - minus);
    indice[idim] = rank / nval;
    rank -= indice[idim] * nval;
  }
}

/**
 * Find the grid node to which the current sample is assigned
 * @param coor   Sample coordinates
 * @param indice Indices of the assigned grid node
 * @param eps    Epsilon to overpass roundoff problem
 * @return
 */
int Grid::coordinateToIndice(const VectorDouble& coor,
                              VectorInt& indice,
                              double eps) const
{
  int ndim = _nDim;
  VectorDouble work1(ndim);
  VectorDouble work2(ndim);

  // Check if all coordinates are defined 

  for (int idim=0; idim<ndim; idim++)
    if (FFFF(coor[idim])) return -1;

  // Shift by the origin

  for (int idim=0; idim<ndim; idim++)
    work1[idim] = coor[idim] - _x0[idim];

  // Perform the Inverse rotation

 _rotation.rotateDirect(work1,work2);

  // Calculate the indices

  for (int idim=0; idim<ndim; idim++)
  {
    int ix = (int) floor(work2[idim] / _dx[idim] + eps);
    if (ix < 0 || ix >= _nx[idim]) return 1;
    indice[idim] = ix;
  }
  return 0;
}

int Grid::coordinateToRank(const VectorDouble& coor, double eps) const
{
  int ndim = _nDim;
  VectorInt indice(ndim);
  if (coordinateToIndice(coor,indice,eps)) return -1;
  return indiceToRank(indice);
}

bool Grid::_isSpaceDimensionValid(int idim) const
{
  if (idim < 0 || idim >= _nDim)
  {
    mesArg("Argument 'idim' is invalid",idim,_nDim);
    return false;
  }
  return true;
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
  int ndim = _nDim;

  // Initiate the Iterator
  _iter = 0;

  // Define the number of positions per space index
  _counts = _nx;

  // Define the order
  if (order.empty() || ndim != (int) order.size())
  {
    _order.resize(ndim,0);
    for (int idim = 0; idim < ndim; idim++)
      _order[idim] = idim;
  }
  else
  {
    // Check that the order array provided is valid
    for (int idim = 0; idim < ndim; idim++)
    {
      bool found = false;
      for (int jdim = 0; jdim < ndim; jdim++)
      {
        int rank = ABS(order[jdim]) - 1;
        if (rank == idim) found = true;
      }
      if (! found)
        my_throw("When provided, 'order' should contain all Space dimensions");
    }
    _order = order;
  }

  // Count the total number of iterations
  _nprod = 1;
  for (int idim = 0; idim < ndim; idim++)
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
  int ndim = _nDim;
  int nval = _nprod;

  VectorInt indices(ndim);
  for (int jdim = ndim - 1; jdim >= 0; jdim--)
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
  int ndim = _nDim;
  VectorInt indg(ndim);
  VectorDouble coor(ndim);

  /* Preliminary checks */

  if (mode != 1 && mode != -1) return;

  /* Get the number of grid nodes */

  for (int idim = 0; idim < ndim; idim++)
  {
    nx[idim] = getNX(idim) + 2 * mode * nshift[idim];
    if (nx[idim] <= 0) return;
    dx[idim] = getDX(idim);
  }

  /* Get the lower left corner of the small grid */

  for (int idim = 0; idim < ndim; idim++)
    indg[idim] = -mode * nshift[idim];
  indicesToCoordinate(indg, coor);

  /* Calculate the center of the lower left cell */

  for (int idim = 0; idim < ndim; idim++)
    x0[idim] = coor[idim];
}

/****************************************************************************/
/*!
 **  Returns the characteristics of a multiple grid
 **
 ** \param[in]  nmult Array of multiplicity coefficients
 ** \param[in]  flag_cell 1 for cell matching; 0 for point matching
 **
 ** \param[out] nx    Array of number of grid meshes
 ** \param[out] dx    Array of grid meshes
 ** \param[out] x0    Array of grid origins
 **
 *****************************************************************************/
void Grid::multiple(const VectorInt& nmult,
                     int flag_cell,
                     VectorInt& nx,
                     VectorDouble& dx,
                     VectorDouble& x0) const
{
  int ndim = _nDim;
  VectorInt indg(ndim);
  VectorDouble perc(ndim);
  VectorDouble coor1(ndim);
  VectorDouble coor2(ndim);

  /* Get the number of grid nodes */

  for (int idim = 0; idim < ndim; idim++)
  {
    double value = (double) getNX(idim);
    if (flag_cell)
      nx[idim] = (int) floor(value / (double) nmult[idim]);
    else
      nx[idim] = 1 + (int) floor((value - 1.) / (double) nmult[idim]);
  }

  /* Get the new grid meshes */

  for (int idim = 0; idim < ndim; idim++)
    dx[idim] = getDX(idim) * nmult[idim];

  /* Get the lower left corner of the small grid */

  for (int idim = 0; idim < ndim; idim++)
    indg[idim] = 0;
  for (int idim = 0; idim < ndim; idim++)
    perc[idim] = -0.5;
  indicesToCoordinateInPlace(indg, coor1, perc);
  for (int idim = 0; idim < ndim; idim++)
    perc[idim] = 0.5;
  indicesToCoordinateInPlace(indg, coor2, perc);

  /* Calculate the center of the lower left cell */

  for (int idim = 0; idim < ndim; idim++)
  {
    double delta = (coor2[idim] - coor1[idim]) / 2.;
    if (flag_cell)
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
 ** \param[in]  flag_cell 1 for cell matching; 0 for point matching
 **
 ** \param[out] nx    Array of number of grid meshes
 ** \param[out] dx    Array of grid meshes
 ** \param[out] x0    Array of grid origins
 **
 *****************************************************************************/
void Grid::divider(const VectorInt& nmult,
                    int flag_cell,
                    VectorInt& nx,
                    VectorDouble& dx,
                    VectorDouble& x0) const
{
  int ndim = _nDim;
  VectorInt indg(ndim);
  VectorDouble perc(ndim);
  VectorDouble coor1(ndim);
  VectorDouble coor2(ndim);

  /* Get the number of grid nodes */

  for (int idim = 0; idim < ndim; idim++)
  {
    if (flag_cell)
      nx[idim] = getNX(idim) * nmult[idim];
    else
      nx[idim] = 1 + (getNX(idim) - 1) * nmult[idim];
  }

  /* Get the new grid meshes */

  for (int idim = 0; idim < ndim; idim++)
    dx[idim] = getDX(idim) / ((double) nmult[idim]);

  /* Get the lower left corner of the small grid */

  for (int idim = 0; idim < ndim; idim++)
    indg[idim] = 0;
  for (int idim = 0; idim < ndim; idim++)
    perc[idim] = -0.5;
  indicesToCoordinateInPlace(indg, coor1, perc);
  for (int idim = 0; idim < ndim; idim++)
    perc[idim] = 0.5;
  indicesToCoordinateInPlace(indg, coor2, perc);

  /* Calculate the center of the lower left cell */

  for (int idim = 0; idim < ndim; idim++)
  {
    double delta = (coor2[idim] - coor1[idim]) / 2.;
    if (flag_cell)
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
  int nx = _nx[idim];
  int nmax = nx - 1;
  while (!(ix >= 0 && ix < nx))
  {
    if (ix < 0)
      ix = -ix;
    else if (ix > nmax) ix = 2 * nmax - ix;
  }
  return (ix);
}
