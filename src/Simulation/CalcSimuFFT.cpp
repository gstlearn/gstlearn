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
#include "Simulation/CalcSimuFFT.hpp"
#include "Basic/Law.hpp"
#include "Core/fftn.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/ModelGeneric.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Simulation/SimuFFTParam.hpp"

#include <cmath>

// #define IND(ix, iy, iz) ((iz) + _dims[2] * ((iy) + _dims[1] * (ix)))
#define IND(ix, iy, iz) ((ix) + _dims[0] * ((iy) + _dims[1] * (iz)))

namespace gstlrn
{
  CalcSimuFFT::CalcSimuFFT(Id nbsimu, bool verbose, Id seed)
    : ACalcSimuGaussian(nbsimu, seed, verbose)
    , _param()
    , _nxyz(0)
    , _nx()
    , _shift()
    , _dims()
    , _dim2()
    , _sizes_alloc(0)
    , _cmat()
    , _rnd()
    , _u()
    , _v()
    , _neigh(nullptr)
    , _c00(0.)
  {
    setFlagOnGridOnly(true);
    _neigh = NeighUnique::create();
    setNeigh(_neigh);
  }

  CalcSimuFFT::~CalcSimuFFT()
  {
    delete _neigh;
  }

  bool CalcSimuFFT::_initializeSimulations()
  {
    /* Construction of the Simu_FFT structure and core allocation */

    _allocate();

    /* Preparation of the FFT environment */

    _prepar();

    return true;
  }

  /**
   * @brief Internal function to calculate the covariance value for an increment
   *
   * @param delta Increment vector
   * @return double
   */
  double CalcSimuFFT::_getIncrCovariance(const VectorDouble& delta) const
  {
    auto ndim = _getNDim();
    VectorDouble incr(ndim);
    for (Id i = 0; i < ndim; ++i)
    {
      incr[i] = delta[i];
    }
    return getModelGeneric()->evaluateOneIncr(incr);
  }

  /****************************************************************************/
  /*!
   **  Dimension the ST_FFT structure
   **
   *****************************************************************************/
  void CalcSimuFFT::_allocate()
  {
    _dbgrid = dynamic_cast<DbGrid*>(getDbout());

    _nx.resize(3, 0);
    _shift.resize(3, 0);
    _dims.resize(3, 0);
    _dim2.resize(3, 0);
    _nxyz = 1;
    auto ndim = _getNDim();

    for (Id i = 0; i < 3; i++)
    {
      if (i < ndim)
      {
        _nx[i] = _dbgrid->getNX(i);
        _nxyz *= _nx[i];
      }
    }

    /* Covariance at the origin (calculated once for all) */

    VectorDouble d0(ndim, 0.);
    _c00 = _getIncrCovariance(d0);

    /* Dilate the grid */

    _gridDilate();

    /* Determine the grid extension */

    for (Id i = 0; i < 3; i++)
    {
      if (i < ndim)
      {
        _nx[i] = _dbgrid->getNX(i);
        auto nval = _getNOptimalEven(_shift[i] + _nx[i]);
        _dims[i] = nval;
        _dim2[i] = nval / 2;
      }
      else
      {
        _dims[i] = 1;
        _dim2[i] = 0;
      }
    }

    Id total = 1;
    for (Id idim = 0; idim < ndim; idim++) total *= _dims[idim];
    _sizes_alloc = total;

    if (getVerbose())
    {
      message("Grid parameters after Optimal Dilation :\n");
      if (ndim >= 1) message("- Number of Nodes along X = %d\n", _dims[0]);
      if (ndim >= 2) message("- Number of Nodes along Y = %d\n", _dims[1]);
      if (ndim >= 3) message("- Number of Nodes along Z = %d\n", _dims[2]);
      message("- Total count of grid nodes = %d\n", _sizes_alloc);
    }

    /* Core allocation */

    _cmat.resize(_sizes_alloc, 0.);
    _rnd.resize(_sizes_alloc, 0.);
    _u.resize(_sizes_alloc, 0.);
    _v.resize(_sizes_alloc, 0.);
  }

  /****************************************************************************/
  /*!
   **  Returns the closest value, larger than the argument, which is
   **  factorized as the product of low factors
   **
   ** \return  Returned number
   **
   ** \param[in]  number input number
   ** \param[in]  largeFactor Maximum value for a Factor
   **
   *****************************************************************************/
  Id CalcSimuFFT::_getNOptimalEven(Id number, Id largeFactor)
  {
    Id local = number;
    if ((local % 2) == 1) local++;

    bool answer = true;
    while (answer)
    {
      VectorInt factors = _getFactors(local);
      Id nfact = static_cast<Id>(factors.size());
      answer = false;
      for (Id i = 0; i < nfact; i++)
        if (factors[i] > largeFactor) answer = 1;
      if (answer) local += 2;
    }
    return (local);
  }

  /****************************************************************************/
  /*!
   **  Get the factor decomposition of a number
   **
   ** \return  Count of active factors
   **
   ** \param[in]  number  number to be decomposed
   **
   *****************************************************************************/
  VectorInt CalcSimuFFT::_getFactors(Id number)
  {
    VectorInt factors;
    Id local = number;
    Id nfact = 0;

    /* Decomposition in multiples of 2 */

    Id j = 2;
    while ((local % j) == 0)
    {
      factors.push_back(j);
      nfact++;
      local /= j;
    }

    /* Decomposition in higher level multiples */

    j = 3;
    do
    {
      while ((local % j) == 0)
      {
        factors.push_back(j);
        nfact++;
        local /= j;
      }
      j += 2;
    } while (j <= local);

    if (nfact <= 0) factors.push_back(1);

    return factors;
  }

  /****************************************************************************/
  /*!
   **  Calculates the grid extension in a given grid direction
   **
   *****************************************************************************/
  void CalcSimuFFT::_gridDilate()
  {
    /* Origin of the grid */

    auto ndim = _getNDim();

    VectorInt indg(3);
    VectorDouble xyz0(3);
    VectorVectorDouble xyz(3);
    for (Id i = 0; i < 3; i++) xyz[i].resize(3);

    /* Location of the elementary end points */

    indg.fill(0.);
    _dbgrid->rankToCoordinatesInPlace(_dbgrid->indiceToRank(indg), xyz0);

    for (Id i = 0; i < 3; i++)
    {
      indg.fill(0.);
      indg[i] = 1;
      _dbgrid->rankToCoordinatesInPlace(_dbgrid->indiceToRank(indg), xyz[i]);
    }

    /* Coordinates of the grid vector in the rotated space */

    for (Id j = 0; j < 3; j++)
      for (Id i = 0; i < 3; i++)
      {
        xyz[j][i] -= xyz0[i];
        if (j >= ndim) xyz[j][i] = 0.;
      }

    /* Evaluate the count of elementary grid mesh (in each direction) */
    /* for the covariance to become negligeable */

    Id ndx = 1;
    Id ndy = 1;
    Id ndz = 1;

    /* Extension along X */

    if (ndim >= 1)
    {
      bool ok_dir = false;
      while (!ok_dir)
      {
        ok_dir = _checkCorrect(xyz, ndx, 0, 0);
        if (!ok_dir) ndx++;
      }
    }

    /* Extension along Y */

    if (ndim >= 2)
    {
      bool ok_dir = false;
      while (!ok_dir)
      {
        ok_dir = _checkCorrect(xyz, 0, ndy, 0);
        if (!ok_dir) ndy++;
      }
    }

    /* Extension along Z */

    if (ndim >= 3)
    {
      bool ok_dir = false;
      while (!ok_dir)
      {
        ok_dir = _checkCorrect(xyz, 0, 0, ndz);
        if (!ok_dir) ndz++;
      }
    }

    /* Ultimate corrections */

    if (ndim < 3) ndz = 0;
    if (ndim < 2) ndy = 0;
    if (ndim < 1) ndx = 0;

    /* Storing arguments */

    _shift[0] = ndx;
    _shift[1] = ndy;
    _shift[2] = ndz;

    /* Optional printout */

    if (getVerbose())
    {
      message("Grid Dilation parameters :\n");
      if (ndim >= 1) message("- Number of Nodes along X = %d\n", ndx);
      if (ndim >= 2) message("- Number of Nodes along Y = %d\n", ndy);
      if (ndim >= 3) message("- Number of Nodes along Z = %d\n", ndz);
    }
  }

  /****************************************************************************/
  /*!
   **  Checks if the covariance is below threshold for tested distance
   **
   ** \return True if the covariance for the current increment is below threshold; 0 otherwise
   **
   ** \param[in]   xyz     Grid increment
   ** \param[in]   ix      Grid index along X
   ** \param[in]   iy      Grid index along Y
   ** \param[in]   iz      Grid index along Z
   **
   *****************************************************************************/
  bool CalcSimuFFT::_checkCorrect(
    const VectorVectorDouble& xyz,
    Id ix,
    Id iy,
    Id iz)
  {
    auto ndim = _getNDim();
    VectorDouble delta(ndim, 0.);
    for (Id i = 0; i < ndim; i++)
      delta[i] = ix * xyz[i][0] + iy * xyz[i][1] + iz * xyz[i][2];
    double value = _getIncrCovariance(delta);

    double percent = MIN(MAX(_param.getPercent(), 1.0), 100.0);
    return (value <= _c00 * percent / 100.);
  }

  /****************************************************************************/
  /*!
   **  Prepares the simulation on a grid using Discrete FFT
   **
   *****************************************************************************/
  void CalcSimuFFT::_prepar()
  {
    VectorDouble del(3, 0.);
    VectorDouble xyz0(3, 0.);
    VectorDouble xyz(3, 0.);
    VectorDouble delta(3, 0.);
    VectorInt indg(3, 0);
    VectorInt jnd(3);
    VectorVectorDouble xyz1(3);
    for (Id i = 0; i < 3; i++) xyz1[i].resize(3, 0.);
    _dbgrid->rankToCoordinatesInPlace(_dbgrid->indiceToRank(indg), xyz0);

    /* Initializations */

    auto ndim = _getNDim();
    VectorDouble cplx(_sizes_alloc, 0.);
    VectorDouble cply(_sizes_alloc, 0.);

    for (Id i = 0; i < 3; i++)
    {
      indg.fill(0);
      indg[i] = 1;
      if (i < ndim)
      {
        _dbgrid->rankToCoordinatesInPlace(_dbgrid->indiceToRank(indg), xyz1[i]);
        for (Id j = 0; j < 3; j++) xyz1[i][j] -= xyz0[j];
        delta[i] = _dbgrid->getDX(i) * _dbgrid->getNX(i);
      }
      else
      {
        for (Id j = 0; j < 3; j++) xyz1[i][j] = 0.;
        delta[i] = 0.;
      }
    }

    // Filling the covariance matrix on the dilated grid

    Id ecr = 0;
    for (Id iz = 0; iz < _dims[2]; iz++)
      for (Id iy = 0; iy < _dims[1]; iy++)
        for (Id ix = 0; ix < _dims[0]; ix++)
        {
          jnd[0] = (ix <= _dim2[0]) ? ix : ix - _dims[0];
          jnd[1] = (iy <= _dim2[1]) ? iy : iy - _dims[1];
          jnd[2] = (iz <= _dim2[2]) ? iz : iz - _dims[2];

          for (Id i = 0; i < 3; i++)
          {
            double proj = 0.;
            for (Id j = 0; j < 3; j++) proj += jnd[j] * xyz1[j][i];
            del[i] = (i < ndim) ? proj : 0.;
          }
          cplx[ecr] = _getIncrCovariance(del);
          cply[ecr] = 0.;
          ecr++;
        }

    /* Perform the Fast Fourier Transform */

    (void)fftn(ndim, _dims.data(), cplx.data(), cply.data(), -1, 1.);

    /* Looking for negative terms */

    double total_plus = 0.;
    double total_moins = 0.;
    for (Id i = 0; i < _sizes_alloc; i++)
    {
      cplx[i] /= static_cast<double>(_sizes_alloc);
      if (cplx[i] < 0)
        total_moins -= cplx[i];
      else
        total_plus += cplx[i];
    }

    /* Correcting positive terms of the spectrum */

    double correc = (total_plus - total_moins) / total_plus;
    if (total_moins > 0)
    {
      for (Id i = 0; i < _sizes_alloc; i++)
        if (cplx[i] < 0.)
          cplx[i] = 0.;
        else
          cplx[i] *= correc;
    }

    /* Converting into amplitude */

    for (Id i = 0; i < _sizes_alloc; i++) cplx[i] = sqrt(cplx[i] / 2.);

    /* Printout statistics */

    if (getVerbose())
    {
      message("Statistics on the Discrete Periodic Covariance\n");
      message("- Total of positive frequencies = %lf\n", total_plus);
      message("- Total of negative frequencies = %lf\n", total_moins);
      message("\n");
    }

    // Returned arguments

    _cmat = cplx;
    _rnd = cply;
  }

  /****************************************************************************/
  /*!
   **  Initiate a vector of random normal values
   **
   *****************************************************************************/
  bool CalcSimuFFT::_simulate(Id isimu)
  {
    DECLARE_UNUSED(isimu);

    for (Id i = 0; i < _sizes_alloc; i++) _u[i] = _cmat[i] * law_gaussian();
    for (Id i = 0; i < _sizes_alloc; i++) _v[i] = _cmat[i] * law_gaussian();

    switch (_getNDim())
    {
      case 1:
        for (Id ix = 0; ix < _dims[0]; ix += _dim2[0]) _setVariance(ix, 0, 0);
        break;

      case 2:
        for (Id iy = 0; iy < _dims[1]; iy += _dim2[1])
          for (Id ix = 0; ix < _dims[0]; ix += _dim2[0])
            _setVariance(ix, iy, 0);
        break;

      case 3:
        for (Id iz = 0; iz < _dims[2]; iz += _dim2[2])
          for (Id iy = 0; iy < _dims[1]; iy += _dim2[1])
            for (Id ix = 0; ix < _dims[0]; ix += _dim2[0])
              _setVariance(ix, iy, iz);
        break;

      default: break;
    }

    _defineSymmetry();

    return true;
  }

  /****************************************************************************/
  /*!
   **  Correct the variance of the spectrum for real U
   **
   ** \param[in]  ix    Cell location along X
   ** \param[in]  iy    Cell location along Y
   ** \param[in]  iz    Cell location along Z
   **
   *****************************************************************************/
  void CalcSimuFFT::_setVariance(Id ix, Id iy, Id iz)
  {
    Id ind = IND(ix, iy, iz);
    _u[ind] *= sqrt(2.0);
    _v[ind] = 0.;
  }

  /****************************************************************************/
  /*!
   **  Operate the symmetry
   **
   *****************************************************************************/
  void CalcSimuFFT::_defineSymmetry(void)
  {

    /* Dispatch according to the space dimension */

    switch (_getNDim())
    {
      case 1: _defineSym1(); break;

      case 2: _defineSym2(0); break;

      case 3: _defineSym3(); break;

      default: break;
    }
  }

  /****************************************************************************/
  /*!
   **  Operate the symmetry for a 1-D space
   **
   *****************************************************************************/
  void CalcSimuFFT::_defineSym1()

  {
    // A(1) and A(N1/2+1) are real
    for (Id ix = 0; ix < _dims[0]; ix += _dim2[0]) _setZero(ix, 0, 0);

    // A(j) = A*(N1-j+2) for j in [2,N1/2]
    for (Id ix = 1; ix < _dim2[0]; ix++)
    {
      Id jx = _dims[0] - ix;
      _setConjugate(ix, 0, 0, jx, 0, 0);
    }
  }

  /****************************************************************************/
  /*!
   **  Operate the symmetry for a 2-D space
   **
   ** \param[in]  iz0   fixed third index
   **
   *****************************************************************************/
  void CalcSimuFFT::_defineSym2(Id iz0)
  {
    // A(1,1), A(N1/2,1), A(1,N2/2) and A(N1/2,N2/2) real
    for (Id iy = 0; iy < _dims[1]; iy += _dim2[1])
      for (Id ix = 0; ix < _dims[0]; ix += _dim2[0]) _setZero(ix, iy, iz0);

    // A(1,k)      = A*(1,N2-k+2)      for k in [2,N2/2]
    // A(N2/2+1,k) = A*(N2/2+1,N2-k+2) for k in [2,N2/2]
    for (Id iy = 1; iy < _dim2[1]; iy++)
      for (Id ix = 0; ix < _dims[0]; ix += _dim2[0])
      {
        Id jy = _dims[1] - iy;
        _setConjugate(ix, iy, iz0, ix, jy, iz0);
      }

    // A(j,1)      = A*(N1-j+2,1)      for j in [2,N1/2]
    // A(j,N2/2+1) = A*(N1-j+2,N2/2+1) for j in [2,N1/2]
    for (Id iy = 0; iy < _dims[1]; iy += _dim2[1])
      for (Id ix = 1; ix < _dim2[0]; ix++)
      {
        Id jx = _dims[0] - ix;
        _setConjugate(ix, iy, iz0, jx, iy, iz0);
      }

    // A(j,k) = A*(N1-j+2,N2-k+2) for j in [2,N1/2] and k in [2,N2/2]
    for (Id iy = 1; iy < _dim2[1]; iy++)
      for (Id ix = 1; ix < _dim2[0]; ix++)
      {
        Id jx = _dims[0] - ix;
        Id jy = _dims[1] - iy;
        _setConjugate(ix, iy, iz0, jx, jy, iz0);
      }

    // A(j,N2-k+2) = A*(N1-j+2,k) for j in [2,N1/2] and k in [2,N2/2]
    for (Id iy = 1; iy < _dim2[1]; iy++)
      for (Id ix = 1; ix < _dim2[0]; ix++)
      {
        Id jx = _dims[0] - ix;
        Id jy = _dims[1] - iy;
        _setConjugate(ix, jy, iz0, jx, iy, iz0);
      }
  }

  /****************************************************************************/
  /*!
   **  Operate the symmetry for a 3-D space
   **
   *****************************************************************************/
  void CalcSimuFFT::_defineSym3()

  {
    // For l=1 or N3/2+1, use the 2-D symmetry
    for (Id iz = 0; iz < _dims[2]; iz += _dim2[2]) _defineSym2(iz);

    // For the other planes:

    // A(1,1,l)           = A*(1,1,N3-l+2)           for l in [2,N3/2]
    // A(N1/2+1,1,l)      = A*(N1/2+1,1,N3-l+2)      for l in [2,N3/2]
    // A(1,N2/2+1,l)      = A*(1,N2/2+1,N3-l+2)      for l in [2,N3/2]
    // A(N1/2+1,N2/2+1,l) = A*(N1/2+1,N2/2+1,N3-l+2) for l in [2,N3/2]
    for (Id iz = 1; iz < _dim2[2]; iz++)
      for (Id iy = 0; iy < _dims[1]; iy += _dim2[1])
        for (Id ix = 0; ix < _dims[0]; ix += _dim2[0])
        {
          Id jz = _dims[2] - iz;
          _setConjugate(ix, iy, iz, ix, iy, jz);
        }

    // A(1,k,l)      = A*(1,N2-k+2,N3-l+2)      for k in [2,N2/2] and l in [2,N3/2]
    // A(N1/2+1,k,l) = A*(N1/2+1,N2-k+2,N3-l+2) for k in [2,N2/2] and l in [2,N3/2]
    for (Id iz = 1; iz < _dim2[2]; iz++)
      for (Id iy = 1; iy < _dim2[1]; iy++)
        for (Id ix = 0; ix < _dims[0]; ix += _dim2[0])
        {
          Id jy = _dims[1] - iy;
          Id jz = _dims[2] - iz;
          _setConjugate(ix, iy, iz, ix, jy, jz);
        }

    // A(1,N2-k+2,l)      = A*(1,k,N3/-l+2)      for k in [2,N2/2] and l in [2,N3/2]
    // A(N1/2+1,N2-k+2,l) = A*(N1/2+1,k,N3/-l+2) for k in [2,N2/2] and l in [2,N3/2]
    for (Id iz = 1; iz < _dim2[2]; iz++)
      for (Id iy = 1; iy < _dim2[1]; iy++)
        for (Id ix = 0; ix < _dims[0]; ix += _dim2[0])
        {
          Id jy = _dims[1] - iy;
          Id jz = _dims[2] - iz;
          _setConjugate(ix, jy, iz, ix, iy, jz);
        }

    // A(j,1,l)      = A*(N1-j+2,1,N3-l+2)      for j in [2,N1/2] and l in [2,N3/2]
    // A(j,N2/2+1,l) = A*(N1-j+2,N2/2+1,N3-l+2) for j in [2,N1/2] and l in [2,N3/2]
    for (Id iz = 1; iz < _dim2[2]; iz++)
      for (Id iy = 0; iy < _dims[1]; iy += _dim2[1])
        for (Id ix = 1; ix < _dim2[0]; ix++)
        {
          Id jx = _dims[0] - ix;
          Id jz = _dims[2] - iz;
          _setConjugate(ix, iy, iz, jx, iy, jz);
        }

    // A(N1-j+2,1,l)      = A*(j,1,N3-l+2)      for j in [2,N1/2] and l in [2,N3/2]
    // A(N1-j+2,N2/2+1,l) = A*(j,N2/2+1,N3-l+2) for j in [2,N1/2] and l in [2,N3/2]
    for (Id iz = 1; iz < _dim2[2]; iz++)
      for (Id iy = 0; iy < _dims[1]; iy += _dim2[1])
        for (Id ix = 1; ix < _dim2[0]; ix++)
        {
          Id jx = _dims[0] - ix;
          Id jz = _dims[2] - iz;
          _setConjugate(jx, iy, iz, ix, iy, jz);
        }

    // A(j,k,l) = A*(N1-j+2,N2-k+2,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
    for (Id iz = 1; iz < _dim2[2]; iz++)
      for (Id iy = 1; iy < _dim2[1]; iy++)
        for (Id ix = 1; ix < _dim2[0]; ix++)
        {
          Id jx = _dims[0] - ix;
          Id jy = _dims[1] - iy;
          Id jz = _dims[2] - iz;
          _setConjugate(ix, iy, iz, jx, jy, jz);
        }

    // A(N1-j+2,N2-k+2,l) = A*(j,k,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
    for (Id iz = 1; iz < _dim2[2]; iz++)
      for (Id iy = 1; iy < _dim2[1]; iy++)
        for (Id ix = 1; ix < _dim2[0]; ix++)
        {
          Id jx = _dims[0] - ix;
          Id jy = _dims[1] - iy;
          Id jz = _dims[2] - iz;
          _setConjugate(jx, jy, iz, ix, iy, jz);
        }

    // A(j,N2-k+2,l) = A*(N2-j+2,k,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
    for (Id iz = 1; iz < _dim2[2]; iz++)
      for (Id iy = 1; iy < _dim2[1]; iy++)
        for (Id ix = 1; ix < _dim2[0]; ix++)
        {
          Id jx = _dims[0] - ix;
          Id jy = _dims[1] - iy;
          Id jz = _dims[2] - iz;
          _setConjugate(ix, jy, iz, jx, iy, jz);
        }

    // A(N1-j+2,k,l) = A*(j,k,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
    for (Id iz = 1; iz < _dim2[2]; iz++)
      for (Id iy = 1; iy < _dim2[1]; iy++)
        for (Id ix = 1; ix < _dim2[0]; ix++)
        {
          Id jx = _dims[0] - ix;
          Id jy = _dims[1] - iy;
          Id jz = _dims[2] - iz;
          _setConjugate(jx, iy, iz, ix, jy, jz);
        }
  }

  /****************************************************************************/
  /*!
   **  Set the imaginary part of a cell to zero
   **
   ** \param[in]  ix    Cell location along X
   ** \param[in]  iy    Cell location along Y
   ** \param[in]  iz    Cell location along Z
   **
   *****************************************************************************/
  void CalcSimuFFT::_setZero(Id ix, Id iy, Id iz)
  {
    Id ind = IND(ix, iy, iz);
    _v[ind] = 0.;
  }

  /****************************************************************************/
  /*!
   **  Set the target cell as the conjugate of the input cell
   **
   ** \param[in]  ix    Input cell location along X
   ** \param[in]  iy    Input cell location along Y
   ** \param[in]  iz    Input cell location along Z
   ** \param[in]  jx    Target cell location along X
   ** \param[in]  jy    Target cell location along Y
   ** \param[in]  jz    Target cell location along Z
   **
   *****************************************************************************/
  void CalcSimuFFT::_setConjugate(Id ix, Id iy, Id iz, Id jx, Id jy, Id jz)
  {
    Id ind1 = IND(ix, iy, iz);
    Id ind2 = IND(jx, jy, jz);
    _u[ind2] = _u[ind1];
    _v[ind2] = -_v[ind1];
  }

  /****************************************************************************/
  /*!
   **  Perform a non-conditional simulation on the grid
   **
   ** \param[in]  db    Db Grid structure
   ** \param[in]  activeArray    Array of active nodes for writing the simulation
   ** \param[in]  tab   Table of parameters for writing the simulation
   **
   *****************************************************************************/
  void CalcSimuFFT::_compute(
    Db* db,
    const VectorBool& activeArray,
    VectorVectorDouble& tab)
  {
    DECLARE_UNUSED(db);

    (void)fftn(_getNDim(), _dims.data(), _u.data(), _v.data(), 1, 1.);

    Id nx = MAX(_nx[0], 1);
    Id ny = MAX(_nx[1], 1);
    Id nz = MAX(_nx[2], 1);

    /* Retrieving the simulation */

    Id ecr = 0;
    for (Id iz = 0; iz < nz; iz++)
      for (Id iy = 0; iy < ny; iy++)
        for (Id ix = 0; ix < nx; ix++, ecr++)
        {
          if (!activeArray[ecr]) continue;
          Id jx = ix + _shift[0];
          Id jy = iy + _shift[1];
          Id jz = iz + _shift[2];
          tab[0][ecr] = _u[IND(jx, jy, jz)];
        }
  }

  /****************************************************************************/
  /*!
   **  Calculate the mean lognormal covariance over the block
   **
   ** \return  Mean lognormal covariance
   **
   ** \param[in]  sigma  Logarithmic variance value
   **
   *****************************************************************************/
  double CalcSimuFFT::_support(double sigma)
  {
    double value = 0.;
    if (isZero(sigma)) return (TEST);

    switch (_getNDim())
    {
      case 1: value = _support1(sigma); break;

      case 2: value = _support2(sigma); break;

      case 3: value = _support3(sigma); break;

      default: break;
    }

    /* Calculate the scale */

    double scale = 1.;
    for (Id idim = 0; idim < _getNDim(); idim++)
      scale *= (_nx[idim] * _nx[idim]);
    value /= scale;

    /* Back-transform into change of support coefficient */

    if (!FFFF(sigma)) value = log(value) / (sigma * sigma);

    return (sqrt(value));
  }

  /****************************************************************************/
  /*!
   **  Calculate the mean lognormal covariance over the block (1-D)
   **
   ** \return  Mean lognormal covariance
   **
   ** \param[in]  sigma  Logarithmic variance value
   **
   *****************************************************************************/
  double CalcSimuFFT::_support1(double sigma)
  {
    double value = 0.;
    for (Id ix = -_nx[0]; ix <= _nx[0]; ix++)
    {
      Id iix = (ix < 0) ? _dims[0] + ix : ix;
      double rho = _rhoSigma(sigma, iix, 0, 0);
      value += (_nx[0] - ABS(ix)) * rho;
    }
    return (value);
  }

  /****************************************************************************/
  /*!
   **  Calculate the mean lognormal covariance over the block (2-D)
   **
   ** \return  Mean lognormal covariance
   **
   ** \param[in]  sigma  Logarithmic variance value
   **
   *****************************************************************************/
  double CalcSimuFFT::_support2(double sigma)
  {
    double value = 0.;
    for (Id ix = -_nx[0]; ix <= _nx[0]; ix++)
      for (Id iy = -_nx[1]; iy <= _nx[1]; iy++)
      {
        Id iix = (ix < 0) ? _dims[0] + ix : ix;
        Id iiy = (iy < 0) ? _dims[1] + iy : iy;
        double rho = _rhoSigma(sigma, iix, iiy, 0);
        value += ((_nx[0] - ABS(ix)) * (_nx[1] - ABS(iy)) * rho);
      }
    return (value);
  }

  /****************************************************************************/
  /*!
   **  Calculate the mean lognormal covariance over the block (3-D)
   **
   ** \return  Mean lognormal covariance
   **
   ** \param[in]  sigma  Logarithmic variance value
   **
   *****************************************************************************/
  double CalcSimuFFT::_support3(double sigma)
  {
    double value = 0.;
    for (Id ix = -_nx[0]; ix <= _nx[0]; ix++)
      for (Id iy = -_nx[1]; iy <= _nx[1]; iy++)
        for (Id iz = -_nx[2]; iz <= _nx[2]; iz++)
        {
          Id iix = (ix < 0) ? _dims[0] + ix : ix;
          Id iiy = (iy < 0) ? _dims[1] + iy : iy;
          Id iiz = (iz < 0) ? _dims[2] + iz : iz;
          double rho = _rhoSigma(sigma, iix, iiy, iiz);
          value +=
            ((_nx[0] - ABS(ix)) * (_nx[1] - ABS(iy)) * (_nx[2] - ABS(iz))
             * rho);
        }
    return (value);
  }

  /****************************************************************************/
  /*!
   **  Calculate the exponential of the scaled correlation
   **
   ** \return  Value of the transformed correlation
   **
   ** \param[in]  sigma  Logarithmic variance value
   ** \param[in]  ix     Index for the discretized covariance along X
   ** \param[in]  iy     Index for the discretized covariance along Y
   ** \param[in]  iz     Index for the discretized covariance along Z
   **
   *****************************************************************************/
  double CalcSimuFFT::_rhoSigma(double sigma, Id ix, Id iy, Id iz)
  {
    double rho = _cmat[IND(ix, iy, iz)];
    if (!FFFF(sigma)) rho = exp(sigma * sigma * rho);
    return (rho);
  }

  bool CalcSimuFFT::_check()
  {
    if (!ACalcSimuGaussian::_check()) return false;

    auto ndim = _getNDim();
    auto nvar = getNVar();
    if (ndim < 1 || ndim > 3)
    {
      messerr("The FFT Method is not a relevant simulation model");
      messerr("for this Space Dimension (%d)", ndim);
      return false;
    }
    if (!getDbout()->isGrid())
    {
      messerr("The argument 'dbout' should be a grid");
      return false;
    }
    if (nvar != 1)
    {
      messerr(" The FFT method is limited to the Monovariate case");
      return false;
    }
    return true;
  }

  bool CalcSimuFFT::_preprocess()
  {
    return ACalcSimuGaussian::_preprocess();
  }

  void CalcSimuFFT::_rollback()
  {
    _cleanVariableDb(1);
  }

  VectorDouble CalcSimuFFT::changeSupport(const VectorDouble& sigma)
  {
    // Allocation

    _allocate();

    /* Preparation of the FFT environment */

    _prepar();

    /* Calculate the correlation matrix (possibly rescaled) */

    (void)fftn(_getNDim(), _dims.data(), _cmat.data(), _rnd.data(), 1, 1.);

    /* Loop on the different lognormal variances */

    Id nval = static_cast<Id>(sigma.size());
    VectorDouble r2val;
    if (nval > 0)
    {
      r2val.resize(nval, TEST);
      for (Id ival = 0; ival < nval; ival++)
        r2val[ival] = _support(sigma[ival]);
    }
    else
    {
      r2val.resize(1, TEST);
      r2val[0] = _support(TEST);
    }
    return r2val;
  }

} // namespace gstlrn
