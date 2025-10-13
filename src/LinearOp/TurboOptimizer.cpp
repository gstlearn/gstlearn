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
#include "LinearOp/TurboOptimizer.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#define GETADR(nrows, irow, icol) ((icol) * (nrows) + (irow))
#define U(i, j)                   (uu[GETADR(TO_ncorner, i, j)])
#define V(i, j)                   (vv[GETADR(TO_ncorner, i, j)])
#define H(i, j)                   (hh[GETADR(TO_ndim, i, j)])
#define W(i, j)                   (ww[GETADR(TO_ndim, i, j)])

namespace gstlrn
{
TurboOptimizer::TurboOptimizer(Id nx,
                               Id ny,
                               double dx,
                               double dy,
                               double x0,
                               double y0,
                               double scale,
                               double sill,
                               Id param,
                               Id flagOne)
  : _isCalculated(false)
  , _nx(nx)
  , _ny(ny)
  , _dx(dx)
  , _dy(dy)
  , _x0(x0)
  , _y0(y0)
  , _scale(scale)
  , _sill(sill)
  , _param(param)
  , _poncif(0)
  , _center(0)
  , _nxred(0)
  , _half(0)
  , _flagOne(flagOne)
  , _Blin()
  , _TildeC_T()
  , _S_T()
  , _Q_T()
{
}

TurboOptimizer::TurboOptimizer(const TurboOptimizer& tbo)
  : _isCalculated(tbo._isCalculated)
  , _nx(tbo._nx)
  , _ny(tbo._ny)
  , _dx(tbo._dx)
  , _dy(tbo._dy)
  , _x0(tbo._x0)
  , _y0(tbo._y0)
  , _scale(tbo._scale)
  , _sill(tbo._sill)
  , _param(tbo._param)
  , _poncif(tbo._poncif)
  , _center(tbo._center)
  , _nxred(tbo._nxred)
  , _half(tbo._half)
  , _flagOne(tbo._flagOne)
  , _TildeC_T(tbo._TildeC_T)
  , _Lambda_T(tbo._Lambda_T)
  , _S_T(tbo._S_T)
  , _Q_T(tbo._Q_T)
{
}

TurboOptimizer& TurboOptimizer::operator=(const TurboOptimizer& tbo)
{
  if (this != &tbo)
  {
    _isCalculated = tbo._isCalculated;
    _nx           = tbo._nx;
    _ny           = tbo._ny;
    _dx           = tbo._dx;
    _dy           = tbo._dy;
    _x0           = tbo._x0;
    _y0           = tbo._y0;
    _scale        = tbo._scale;
    _sill         = tbo._sill;
    _param        = tbo._param;
    _poncif       = tbo._poncif;
    _center       = tbo._center;
    _nxred        = tbo._nxred;
    _half         = tbo._half;
    _flagOne      = tbo._flagOne;
    _TildeC_T     = tbo._TildeC_T;
    _Lambda_T     = tbo._Lambda_T;
    _S_T          = tbo._S_T;
    _Q_T          = tbo._Q_T;
  }
  return *this;
}

TurboOptimizer::~TurboOptimizer()
{
}

/**
 * Definition of the 2-D non rotated Grid
 * @param nx Number of nodes along X
 * @param ny Number of nodes along Y
 * @param dx Mesh of the grid along X
 * @param dy Mesh of the grid along Y
 * @param x0 Origin of the grid along X
 * @param y0 Origin of the grid along Y
 */
void TurboOptimizer::setGrid(Id nx,
                             Id ny,
                             double dx,
                             double dy,
                             double x0,
                             double y0)
{
  _nx = nx;
  _ny = ny;
  _dx = dx;
  _dy = dy;
  _x0 = x0;
  _y0 = y0;
}

/**
 * Definition of the Model (single isotropic Matérn structure) by range
 * @param range Range of the structure
 * @param sill  Sill of the structure
 * @param param Matérn parameter (third parameter)
 */
void TurboOptimizer::setModelByRange(double range, double sill, Id param)
{
  _sill  = sill;
  _param = param;
  _scale = _rangeToScale(range); // Performed after setting 'param'
}

/**
 * Definition of the Model (single isotropic Matérn structure) by scale
 * @param scale Scale of the structure
 * @param sill  Sill of the structure
 * @param param Matérn parameter (third parameter)
 */
void TurboOptimizer::setModelByScale(double scale, double sill, Id param)
{
  _scale = scale;
  _sill  = sill;
  _param = param;
}

/**
 * Generic method to set all the remaining terms of the class
 * @param flagOne Starting value for numbering of rows and columns in TripletND
 */
void TurboOptimizer::setEnviron(Id flagOne)
{
  _flagOne = flagOne;
}

double TurboOptimizer::_rangeToScale(double range) const
{
  double factor = sqrt(12. * _param);
  return range / factor;
}

Id TurboOptimizer::_getVertex(Id imesh, Id rank) const
{
  VectorInt indice(2);
  Id node, icas;
  _fromMeshToIndex(imesh, &node, &icas);
  _rankToIndice(node, indice, false);
  for (Id idim = 0; idim < TO_ndim; idim++)
    indice[idim] += _MSS(icas, rank, idim);
  return _indiceToRank(indice);
}

void TurboOptimizer::_fromMeshToIndex(Id imesh, Id* node, Id* icas) const
{
  VectorInt indice(2);
  Id rank = imesh / TO_npercell;
  _rankToIndice(rank, indice, true);
  *icas = imesh - rank * TO_npercell;
  *node = _indiceToRank(indice);
}

Id TurboOptimizer::_MSS(Id icas, Id icorn, Id idim0)
{
  static Id S2D[2][3][2] = {{{0, 0}, {1, 0}, {0, 1}}, {{0, 1}, {1, 0}, {1, 1}}};
  return S2D[icas][icorn][idim0];
}

void TurboOptimizer::_rankToIndice(Id rank, VectorInt& indice, bool minusOne) const
{
  if (static_cast<Id>(indice.size()) < 2)
    my_throw("Argument indice should have the correct size");
  Id nval   = minusOne ? (_ny - 1) : _ny;
  indice[1] = rank / nval;
  rank -= indice[1] * nval;
  indice[0] = rank;
}

Id TurboOptimizer::_indiceToRank(VectorInt& indice, bool flag_complete) const
{
  if (flag_complete)
  {
    if (indice[0] < 0 || indice[0] >= _nx) my_throw("Error in indice[0]");
    if (indice[1] < 0 || indice[1] >= _ny) my_throw("Error in indice[1]");
    return _nx * indice[1] + indice[0];
  }
  if (indice[0] < 0 || indice[0] >= _nxred) indice[0] = _center;
  if (indice[1] < 0 || indice[1] >= _nxred) indice[1] = _center;
  return _nxred * indice[1] + indice[0];
}

double TurboOptimizer::_indiceToCoordinate(Id idim0,
                                           const VectorInt& indice) const
{
  if (idim0 == 0) return (_x0 + indice[0] * _dx);
  return (_y0 + indice[1] * _dy);
}

double TurboOptimizer::_getCoor(Id node, Id idim0) const
{
  VectorInt indice(TO_ndim);
  _rankToIndice(node, indice, false);
  return _indiceToCoordinate(idim0, indice);
}

double TurboOptimizer::_getCoorByMesh(Id imesh, Id rank, Id idim0) const
{
  VectorInt indice(TO_ndim);
  auto node = _getVertex(imesh, rank);
  return _getCoor(node, idim0);
}

void TurboOptimizer::_printVector(const std::string& title,
                                  VectorDouble& uu,
                                  I32 width,
                                  I32 ndec)
{
  Id nval = static_cast<Id>(uu.size());
  std::cout << title << std::endl;
  for (Id i = 0; i < nval; i++)
  {
    std::cout << "[" << std::setw(2) << i + 1 << ",] ";
    if (width > 0) std::cout << std::setw(width);
    if (ndec > 0) std::cout << std::setprecision(ndec);
    std::cout << uu[i] << std::endl;
  }
}

void TurboOptimizer::_printMatrix(const std::string& title,
                                  Id nrow,
                                  Id ncol,
                                  VectorDouble& uu,
                                  Id nper_batch,
                                  Id row_shift,
                                  Id col_shift,
                                  I32 width,
                                  I32 ndec)
{
  // Initializations
  Id nbatch = 1 + (ncol - 1) / nper_batch;

  // Printout
  std::cout << title << std::endl;
  for (Id ibatch = 0; ibatch < nbatch; ibatch++)
  {
    Id ncol_min = ibatch * nper_batch;
    Id ncol_max = std::min((ibatch + 1) * nper_batch, ncol);

    std::cout << "     ";
    for (Id icol = ncol_min; icol < ncol_max; icol++)
    {
      std::cout << std::setw(width - 2) << "[," << icol + col_shift + 1 << "]";
    }
    std::cout << std::endl;

    for (Id irow = 0; irow < nrow; irow++)
    {
      std::cout << "[" << std::setw(3) << irow + row_shift + 1 << ",]";
      for (Id icol = ncol_min; icol < ncol_max; icol++)
      {
        std::cout << " ";
        if (width > 0) std::cout << std::setw(width);
        if (ndec > 0) std::cout << std::setprecision(ndec);
        std::cout << uu[GETADR(nrow, irow, icol)];
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void TurboOptimizer::_invert_3x3(VectorDouble& uu,
                                 VectorDouble& vv,
                                 double tol)
{
  double det00 = U(1, 1) * U(2, 2) - U(2, 1) * U(1, 2);
  double det01 = U(1, 0) * U(2, 2) - U(2, 0) * U(1, 2);
  double det02 = U(1, 0) * U(2, 1) - U(2, 0) * U(1, 1);
  double det   = U(0, 0) * det00 - U(0, 1) * det01 + U(0, 2) * det02;
  double toto  = ABS(det);
  if (toto < tol) my_throw("Matrix is not invertible");
  if (ABS(det) < tol) my_throw("Matrix is not invertible");

  V(0, 0) = (U(1, 1) * U(2, 2) - U(2, 1) * U(1, 2)) / det;
  V(0, 1) = -(U(0, 1) * U(2, 2) - U(2, 1) * U(0, 2)) / det;
  V(0, 2) = (U(0, 1) * U(1, 2) - U(1, 1) * U(0, 2)) / det;
  V(1, 0) = -(U(1, 0) * U(2, 2) - U(2, 0) * U(1, 2)) / det;
  V(1, 1) = (U(0, 0) * U(2, 2) - U(2, 0) * U(0, 2)) / det;
  V(1, 2) = -(U(0, 0) * U(1, 2) - U(1, 0) * U(0, 2)) / det;
  V(2, 0) = (U(1, 0) * U(2, 1) - U(2, 0) * U(1, 1)) / det;
  V(2, 1) = -(U(0, 0) * U(2, 1) - U(2, 0) * U(0, 1)) / det;
  V(2, 2) = (U(0, 0) * U(1, 1) - U(1, 0) * U(0, 1)) / det;
}

/**
 * Perform C = A * B (A, B and C are square matrices of dimension 'size')
 * @param size Dimension of the square matrices
 * @param aa First square matrix
 * @param bb Second square matrix
 * @param cc Resulting square matrix
 */
void TurboOptimizer::_prodMatrix(Id size,
                                 const VectorDouble& aa,
                                 const VectorDouble& bb,
                                 VectorDouble& cc)
{
  for (Id i = 0; i < size; i++)
    for (Id j = 0; j < size; j++)
    {
      double value = 0.;
      for (Id k = 0; k < size; k++)
        value += aa[GETADR(size, i, k)] * bb[GETADR(size, k, j)];
      cc[GETADR(size, i, j)] = value;
    }
}

void TurboOptimizer::_prodMatrixVector(Id size,
                                       const VectorDouble& aa,
                                       const VectorDouble& bb,
                                       VectorDouble& cc)
{
  for (Id i = 0; i < size; i++)
  {
    double value = 0.;
    for (Id j = 0; j < size; j++)
      value += aa[GETADR(size, i, j)] * bb[j];
    cc[i] = value;
  }
}

VectorDouble TurboOptimizer::_buildS(const VectorDouble& TildeC) const
{
  Id nvertex_red = _getNVertices_red();
  VectorDouble ss(nvertex_red * nvertex_red, 0.);
  VectorDouble hh(TO_ndim * TO_ndim);
  VectorDouble uu(TO_ncorner * TO_ncorner);
  VectorDouble vv(TO_ncorner * TO_ncorner);
  VectorDouble ww(TO_ndim * TO_ncorner);

  // Load the matrix HH which accounts for the Model

  _loadHH(hh);

  /* Loop on the meshes */

  for (Id imesh = 0; imesh < _getNMeshes(); imesh++)
  {
    for (Id icorn = 0; icorn < TO_ncorner; icorn++)
    {
      // Filling the 3x3 matrix
      for (Id idim = 0; idim < TO_ndim; idim++)
        V(idim, icorn) = _getCoorByMesh(imesh, icorn, idim);
      V(TO_ndim, icorn) = 1.;
    }

    // Inverting matrix
    _invert_3x3(vv, uu);

    // Transposing and compressing the matrix
    for (Id icorn = 0; icorn < TO_ncorner; icorn++)
      for (Id idim = 0; idim < TO_ndim; idim++)
        W(idim, icorn) = U(icorn, idim);

    for (Id irow = 0; irow < TO_ncorner; irow++)
      for (Id icol = 0; icol < TO_ncorner; icol++)
      {
        double value = 0.;
        for (Id k = 0; k < TO_ndim; k++)
          for (Id l = 0; l < TO_ndim; l++)
            value += W(k, irow) * H(k, l) * W(l, icol);

        auto ip1 = _getVertex(imesh, irow);
        auto ip2 = _getVertex(imesh, icol);
        ss[GETADR(nvertex_red, ip1, ip2)] += _getMeshSize() * value;
      }
  }

  // Normalize S matrix by TildeC
  for (Id ip = 0; ip < nvertex_red; ip++)
    for (Id jp = 0; jp < nvertex_red; jp++)
      ss[GETADR(nvertex_red, ip, jp)] /= sqrt(TildeC[ip] * TildeC[jp]);

  return ss;
}

Id TurboOptimizer::_determineInternalGrid(bool verbose)
{
  Id nblin = _param + 2;
  _half    = nblin - 1;
  _poncif  = 2 * _half + 1;
  _nxred   = (1 + _half) + _poncif + (1 + _half);
  _center  = (1 + _half) + _half;

  if (verbose)
  {
    std::cout << "Internal Grid Determination" << std::endl;
    std::cout << "- Matérn parameter = " << _param << std::endl;
    std::cout << "- Dimension of Internal Square Grid = " << _nxred << std::endl;
  }

  if (_nx < _nxred || _ny < _nxred)
  {
    std::cout << "The output grid must be larger than the internal one ("
              << _nxred << ")" << std::endl;
    return 1;
  }
  return 0;
}

/**
 * This function is compulsory as it performs the calculations
 * and allows retrieval of the matrices
 * @param verbose Verbose flag
 */
void TurboOptimizer::run(bool verbose)
{
  if (_isCalculated) return;
  Id nx_memo = _nx;
  Id ny_memo = _ny;

  // Modify the Grid dimension to the Template one
  if (_determineInternalGrid(verbose)) my_throw("Incompatible Grid sizes");
  _nx            = _nxred;
  _ny            = _nxred;
  Id nvertex_red = _getNVertices();

  if (verbose)
    std::cout << "Scale = " << _scale << std::endl;

  _Blin = _buildBlin();
  if (verbose) _printVector("Template Blin Vector", _Blin);

  // Build TildeC template vector
  _TildeC_T = _buildTildeC();
  if (verbose) _printVector("Template TildeC Vector", _TildeC_T);

  // Build Lambda template vector
  _Lambda_T = _buildLambda(_TildeC_T);
  if (verbose) _printVector("Template Lambda Vector", _Lambda_T);

  // Build S template matrix
  _S_T = _buildS(_TildeC_T);
  if (verbose) _printMatrix("Template S Matrix", nvertex_red, nvertex_red, _S_T, 7);

  // Build Q matrix
  _Q_T = _buildQ(_S_T, _Blin, _Lambda_T);
  if (verbose) _printMatrix("Template Q Matrix", nvertex_red, nvertex_red, _Q_T, 6);

  // Restore the initial Grid dimensions
  _nx           = nx_memo;
  _ny           = ny_memo;
  _isCalculated = true;
}

VectorDouble TurboOptimizer::_buildTildeC() const
{
  auto nvertex = _getNVertices();
  VectorDouble TildeC(nvertex, 0.);

  /* Loop on the meshes */

  for (Id imesh = 0; imesh < _getNMeshes(); imesh++)
  {
    for (Id icorn = 0; icorn < TO_ncorner; icorn++)
    {
      auto jp = _getVertex(imesh, icorn);
      TildeC[jp] += _getMeshSize();
    }
  }

  /* Scaling */

  for (Id ip = 0; ip < nvertex; ip++)
    TildeC[ip] /= TO_ncorner;

  return TildeC;
}

VectorDouble TurboOptimizer::_buildLambda(const VectorDouble& TildeC) const
{
  auto nvertex = _getNVertices();
  VectorDouble Lambda(nvertex, 0);
  double value = _scale * _scale;
  for (Id ip = 0; ip < _getNVertices(); ip++)
    Lambda[ip] = sqrt(TildeC[ip] / (value * _sill));

  return Lambda;
}

VectorDouble TurboOptimizer::_buildBlin() const
{
  double PI = 3.14159265358979323846;
  double v1, v2;

  Id ndims2 = static_cast<Id>(TO_ndim / 2.);
  Id p      = _param + ndims2;

  double gammap = 0.;
  for (Id i = 1; i < _param; i++)
    gammap += log(static_cast<double>(i));
  gammap        = exp(gammap);
  double gammaa = 0.;
  for (Id i = 1; i < _param + ndims2; i++)
    gammaa += log(static_cast<double>(i));
  gammaa = exp(gammaa);

  double g0     = pow(4. * PI, ndims2);
  double correc = gammap / (g0 * gammaa);

  VectorDouble blin(p + 1, 0.);
  for (Id i = 0; i <= p; i++)
  {
    // Calculate cnp(p,i)
    v1 = v2 = 0.;
    for (Id j = 0; j < i; j++)
    {
      v1 += log(p - j);
      v2 += log(j + 1);
    }
    double cnp = exp(v1 - v2);
    blin[i]    = cnp * correc;
  }
  return blin;
}

VectorDouble TurboOptimizer::_buildQ(const VectorDouble& ss,
                                     const VectorDouble& blin,
                                     const VectorDouble& lambda) const
{
  auto nvertex = _getNVertices();
  Id nblin     = static_cast<Id>(blin.size());
  VectorDouble qq(nvertex * nvertex, 0.);
  VectorDouble bi(nvertex * nvertex, 0.);
  VectorDouble be(nvertex * nvertex, 0.);

  // First term
  for (Id i = 0; i < nvertex; i++)
    qq[GETADR(nvertex, i, i)] = blin[0];
  for (Id i = 0; i < nvertex * nvertex; i++)
    bi[i] = ss[i];

  // Loop on the different terms
  for (Id iterm = 1; iterm < nblin; iterm++)
  {
    for (Id i = 0; i < nvertex * nvertex; i++)
      qq[i] += bi[i] * blin[iterm];
    if (iterm < nblin - 1)
    {
      _prodMatrix(nvertex, ss, bi, be);
      for (Id i = 0; i < nvertex * nvertex; i++) bi[i] = be[i];
    }
  }

  // Final scaling
  for (Id i = 0; i < nvertex; i++)
    for (Id j = 0; j < nvertex; j++)
      qq[GETADR(nvertex, i, j)] *= lambda[i] * lambda[j];

  return qq;
}

void TurboOptimizer::_loadHH(VectorDouble& hh) const
{
  double value = _scale * _scale;
  H(0, 0)      = value;
  H(0, 1)      = 0.;
  H(1, 0)      = 0.;
  H(1, 1)      = value;
}

/**
 * Display the parameters of the Method (Grid and Model parameters)
 */
void TurboOptimizer::printClass() const
{
  std::cout << "Grid Definition" << std::endl;
  std::cout << "NX = " << _nx << std::endl;
  std::cout << "NY = " << _ny << std::endl;
  std::cout << "DX = " << _dx << std::endl;
  std::cout << "DY = " << _dy << std::endl;
  std::cout << "X0 = " << _x0 << std::endl;
  std::cout << "Y0 = " << _y0 << std::endl;
  std::cout << std::endl;
  std::cout << "Model Definition" << std::endl;
  std::cout << "Scale = " << _scale << std::endl;
  std::cout << "Sill  = " << _sill << std::endl;
  std::cout << "Param = " << _param << std::endl;
  std::cout << "Triplet Numbering starting value = " << _flagOne << std::endl;
}

/**
 * Print the elements of the Internal Meshing
 */
void TurboOptimizer::printMeshes() const
{
  std::cout << "Number of Meshes   '       = " << _getNMeshes() << std::endl;
  std::cout << "Number of Corners per Mesh = " << TO_ncorner << std::endl;
  std::cout << "Number of Vertices         = " << _getNVertices() << std::endl;
  std::cout << "Number of Coordinates      = " << TO_ndim << std::endl;

  for (Id imesh = 0; imesh < _getNMeshes(); imesh++)
  {
    std::cout << "Mesh #" << imesh + 1 << " : ";
    for (Id ic = 0; ic < TO_ncorner; ic++)
      std::cout << _getVertex(imesh, ic) << " ";
    std::cout << std::endl;
  }
  std::cout << std::endl;

  for (Id node = 0; node < _getNVertices(); node++)
  {
    std::cout << "Vertex #" << node + 1 << " : ";
    for (Id idim = 0; idim < TO_ndim; idim++)
      std::cout << _getCoor(node, idim) << " ";
    std::cout << std::endl;
  }
}

void TurboOptimizer::_updateMargin(Id idim0, VectorInt& indice) const
{
  Id nmax = (idim0 == 0) ? _nx : _ny;
  if (indice[idim0] < _half) return;
  if ((nmax - 1) - indice[idim0] < _half)
  {
    indice[idim0] = (_nxred - 1) - ((nmax - 1) - indice[idim0]);
  }
  else
  {
    indice[idim0] = _center;
  }
}

void TurboOptimizer::_getRankInTemplate(VectorInt& indice1,
                                        VectorInt& indice2) const
{
  Id decalx = indice2[0] - indice1[0];
  Id decaly = indice2[1] - indice1[1];
  _updateMargin(0, indice1);
  _updateMargin(1, indice1);
  _updateMargin(0, indice2);
  _updateMargin(1, indice2);
  indice2[0] = indice1[0] + decalx;
  indice2[1] = indice1[1] + decaly;
}

VectorDouble TurboOptimizer::_getVectorFromTemplate(const VectorDouble& vecin) const
{
  auto nvertex = _getNVertices();
  VectorDouble vecout(nvertex, 0.);
  if (!_isCalculated)
    my_throw("You must use the method 'run' beforehand");

  VectorInt indice(TO_ndim, 0);

  for (Id ix = 0; ix < _nx; ix++)
    for (Id iy = 0; iy < _ny; iy++)
    {
      indice[0] = ix;
      indice[1] = iy;
      Id ecr    = _indiceToRank(indice);
      _updateMargin(0, indice);
      _updateMargin(1, indice);
      Id lec      = _indiceToRank(indice, false);
      vecout[ecr] = vecin[lec];
    }
  return vecout;
}

TripletND TurboOptimizer::_getMatrixFromTemplate(const VectorDouble& matin,
                                                 Id nperline) const
{
  VectorInt indice1(TO_ndim, 0);
  VectorInt indice2(TO_ndim, 0);

  // Pre-dimensioning of TripletND structure (to avoid iterative concatenations
  // Trying to guess the number of non-zero terms in matrices

  Id estimated_size = _getNVertices() * nperline;
  TripletND triplet;
  triplet.rows.resize(estimated_size);
  triplet.cols.resize(estimated_size);
  triplet.values.resize(estimated_size);
  if (!_isCalculated)
    my_throw("You must use the method 'run' beforehand");

  // Loop on the vertices (first point)

  Id added_element = 0;
  for (Id iy1 = 0; iy1 < _ny; iy1++)
    for (Id ix1 = 0; ix1 < _nx; ix1++)
    {

      // Loop on the vertices (second point)

      Id ix2min = MAX(0, ix1 - _poncif);
      Id ix2max = MIN(_nx, ix1 + _poncif);
      Id iy2min = MAX(0, iy1 - _poncif);
      Id iy2max = MIN(_ny, iy1 + _poncif);

      for (Id iy2 = iy2min; iy2 < iy2max; iy2++)
        for (Id ix2 = ix2min; ix2 < ix2max; ix2++)
        {

          // Discard comparison between two distant pixels (protection code)
          if (ABS(ix1 - ix2) > _poncif) continue;
          if (ABS(iy1 - iy2) > _poncif) continue;

          indice1[0] = ix1;
          indice1[1] = iy1;
          indice2[0] = ix2;
          indice2[1] = iy2;
          Id ecr1    = _indiceToRank(indice1, true);
          Id ecr2    = _indiceToRank(indice2, true);
          _getRankInTemplate(indice1, indice2);
          Id lec1 = _indiceToRank(indice1, false);
          Id lec2 = _indiceToRank(indice2, false);

          double value = matin[GETADR(_getNVertices_red(), lec1, lec2)];
          if (value != 0.)
          {
            triplet.rows[added_element]   = (ecr1 + _flagOne);
            triplet.cols[added_element]   = (ecr2 + _flagOne);
            triplet.values[added_element] = value;
            added_element++;
            if (added_element >= estimated_size)
              my_throw("Reconsider the pre-estimation of matrix dimensions");
          }
        }
    }

  // Final resizing
  if (estimated_size != added_element)
  {
    triplet.rows.resize(added_element);
    triplet.cols.resize(added_element);
    triplet.values.resize(added_element);
  }
  return triplet;
}

/**
 * Allows retrieving the vector Blin
 * @return The Blin vector
 */
VectorDouble TurboOptimizer::getBlin() const
{
  if (!_isCalculated)
    my_throw("You must use the method 'run' beforehand");
  return _Blin;
}

/**
 * Allows retrieving the vector TildeC
 * @return The TildeC vector
 */
VectorDouble TurboOptimizer::getTildeC() const
{
  return _getVectorFromTemplate(_TildeC_T);
}

/**
 * Allows retrieving the vector Lambda
 * @return The Lambda vector
 */

VectorDouble TurboOptimizer::getLambda() const
{
  return _getVectorFromTemplate(_Lambda_T);
}

/**
 * Allows retrieving the S sparse matrix
 * @return The returned matrix stored as TripletNDs
 */
TripletND TurboOptimizer::getS() const
{
  Id nperline = 5 * _param;
  return _getMatrixFromTemplate(_S_T, nperline);
}

/**
 * Allows retrieving the Q sparse matrix
 * @return The returned matrix stored as Triplets
 */
TripletND TurboOptimizer::getQ() const
{
  Id nbp1     = static_cast<Id>(_Blin.size()) - 1;
  Id nperline = 4 * nbp1 * nbp1;
  return _getMatrixFromTemplate(_Q_T, nperline);
}

Id TurboOptimizer::_coordinateToIndice(double x,
                                       double y,
                                       VectorInt& indice) const
{
  if (static_cast<Id>(indice.size()) < 2)
    my_throw("Argument indice should have the correct size");
  Id ix = static_cast<Id>(floor((x - _x0) / _dx));
  if (ix < 0 || ix >= _nx) return 1;
  indice[0] = ix;
  Id iy     = static_cast<Id>(floor((y - _y0) / _dy));
  if (iy < 0 || iy >= _ny) return 1;
  indice[1] = iy;
  return 0;
}

Id TurboOptimizer::_addWeights(Id icas,
                               double x,
                               double y,
                               const VectorInt& indg0,
                               VectorInt& indices,
                               VectorDouble& lambda) const
{
  VectorDouble lhs(TO_ncorner * TO_ncorner);
  VectorDouble lhsinv(TO_ncorner * TO_ncorner);
  VectorDouble rhs(TO_ncorner);
  VectorInt indgg(2);

  for (Id icorner = 0; icorner < TO_ncorner; icorner++)
  {
    // Generate the indices of the mesh apex
    for (Id idim = 0; idim < TO_ndim; idim++)
      indgg[idim] = indg0[idim] + _MSS(icas, icorner, idim);
    if (indgg[0] < 0 || indgg[0] >= _nx) return 1;
    if (indgg[1] < 0 || indgg[1] >= _ny) return 1;
    indices[icorner] = _indiceToRank(indgg);

    // Update the LHS matrix
    for (Id idim = 0; idim < TO_ndim; idim++)
      lhs[GETADR(TO_ncorner, idim, icorner)] = _indiceToCoordinate(idim, indgg);
    lhs[GETADR(TO_ncorner, TO_ndim, icorner)] = 1.;
  }

  // Generate the right-hand side
  rhs[0] = x;
  rhs[1] = y;
  rhs[2] = 1.;

  // Invert the matrix
  _invert_3x3(lhs, lhsinv);

  // Calculate the weights
  _prodMatrixVector(TO_ncorner, lhsinv, rhs, lambda);

  // Check that all weights are positive
  for (Id idim = 0; idim < TO_ndim; idim++)
    if (lambda[idim] < 0) return 1;

  return 0;
}

/**
 * Returns the weights for interpolating points on the meshing
 * @param x       Vector of X-coordinates for the target points
 * @param y       Vector of Y-coordinates for the target points
 * @return The triplet structure giving the interpolation weights
 * @note: The triplets contain:
 * @note: - rows: the index of the target point
 * @note: - cols: the index of the node of the grid
 * @note: - values: the corresponding weight
 */
TripletND TurboOptimizer::interpolate(const VectorDouble& x,
                                      const VectorDouble& y) const
{
  VectorInt indg0(2);
  VectorInt indices(TO_ncorner);
  VectorDouble lambda(TO_ncorner);
  Id nech = static_cast<Id>(x.size());

  // Pre-allocate the triplet structure

  TripletND triplet;
  Id size = TO_ncorner * nech;
  triplet.rows.resize(size, -1);
  triplet.cols.resize(size, -1);
  triplet.values.resize(size, 0.);

  /* Loop on the samples */

  Id ecr = 0;
  for (Id iech = 0; iech < nech; iech++)
  {

    /* Calculate the grid indices */

    if (_coordinateToIndice(x[iech], y[iech], indg0) == 0)
    {

      /* Loop on the different meshes constituting the cell */

      Id found = -1;
      for (Id icas = 0; icas < TO_npercell && found < 0; icas++)
      {
        if (_addWeights(icas, x[iech], y[iech], indg0, indices, lambda) == 0)
        {
          for (Id icorner = 0; icorner < TO_ncorner; icorner++)
          {
            triplet.rows[ecr + icorner]   = iech + _flagOne;
            triplet.cols[ecr + icorner]   = indices[icorner] + _flagOne;
            triplet.values[ecr + icorner] = lambda[icorner];
          }
          found = icas;
        }
      }
    }
    ecr += TO_ncorner;
  }

  return triplet;
}

/**
 * Expand a sub-part of a Sparse matrix stored as triplets
 * @param row_begin  Starting Row number (included) of the Matrix to be expanded
 * @param row_end    Ending Row number (included) of the Matrix to be expanded
 * @param col_begin  Starting Column number (included) of the Matrix to be expanded
 * @param col_end    Ending Column number (included) of the Matrix to be expanded
 * @param triplet Input matrix stored as TripletNDs
 * @return Matrix stored in full format
 */
VectorDouble TurboOptimizer::_expandTripletToMatrix(Id row_begin,
                                                    Id row_end,
                                                    Id col_begin,
                                                    Id col_end,
                                                    const TripletND& triplet) const
{
  Id nrows = row_end - row_begin + 1;
  Id ncols = col_end - col_begin + 1;
  VectorDouble mat(nrows * ncols, 0.);
  Id size = static_cast<Id>(triplet.rows.size());

  for (Id i = 0; i < size; i++)
  {
    Id irow = triplet.rows[i] - _flagOne;
    Id icol = triplet.cols[i] - _flagOne;
    if (irow >= row_begin && irow <= row_end &&
        icol >= col_begin && icol <= col_end)
      mat[GETADR(nrows, irow - row_begin, icol - col_begin)] = triplet.values[i];
  }
  return mat;
}

void TurboOptimizer::printS(Id nper_batch,
                            Id row_begin,
                            Id row_end,
                            Id col_begin,
                            Id col_end) const
{
  // Convert from IHM to internal numbering
  if (row_begin > 0) row_begin--;
  if (row_end > 0) row_end--;
  if (col_begin > 0) col_begin--;
  if (col_end > 0) col_end--;
  if (row_end <= 0 || row_end >= _getNVertices()) row_end = _getNVertices() - 1;
  if (col_end <= 0 || col_end >= _getNVertices()) col_end = _getNVertices() - 1;

  VectorDouble temp = _expandTripletToMatrix(row_begin, row_end,
                                             col_begin, col_end, getS());
  Id nrows          = row_end - row_begin + 1;
  Id ncols          = col_end - col_begin + 1;
  _printMatrix("Matrix S", nrows, ncols, temp, nper_batch, row_begin, col_begin);
}

void TurboOptimizer::printQ(Id nper_batch,
                            Id row_begin,
                            Id row_end,
                            Id col_begin,
                            Id col_end) const
{
  // Convert from IHM to internal numbering
  if (row_begin > 0) row_begin--;
  if (row_end > 0) row_end--;
  if (col_begin > 0) col_begin--;
  if (col_end > 0) col_end--;
  if (row_end <= 0 || row_end >= _getNVertices()) row_end = _getNVertices() - 1;
  if (col_end <= 0 || col_end >= _getNVertices()) col_end = _getNVertices() - 1;

  VectorDouble temp = _expandTripletToMatrix(row_begin, row_end,
                                             col_begin, col_end, getQ());
  Id nrows          = row_end - row_begin + 1;
  Id ncols          = col_end - col_begin + 1;
  _printMatrix("Matrix Q", nrows, ncols, temp, nper_batch, row_begin, col_begin);
}
} // namespace gstlrn
