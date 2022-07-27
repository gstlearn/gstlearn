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
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Simulation/ACalcSimulation.hpp"
#include "Simulation/SimuFFT.hpp"
#include "Simulation/SimuFFTParam.hpp"
#include "Basic/Law.hpp"

#include <math.h>

#define IND(ix,iy,iz) ((iz) + _dims[2] * ((iy) + _dims[1] * (ix)))
#define U(ix,iy,iz)   (_u[IND(ix,iy,iz)])

SimuFFT::SimuFFT(int nbsimu, int seed)
    : ACalcSimulation(nbsimu, seed),
      _ndim(0),
      _nxyz(0),
      _nx(),
      _shift(),
      _dims(),
      _dim2(),
      _sizes_alloc(0),
      _cmat(),
      _rnd(),
      _u(),
      _v()
{
}

SimuFFT::~SimuFFT()
{
}

int SimuFFT::simulate(DbGrid *db,
                      Model* model,
                      const SimuFFTParam& param,
                      int iptr,
                      bool verbose)
{
  /* Define the global variables */

  law_set_random_seed(getSeed());

  // Check validity of the context

  if (! _isValid(db, model)) return 1;

  /* Construction of the Simu_FFT structure and core allocation */

  _alloc(db, model, param, verbose);

  /* Preparation of the FFT environment */

  _prepar(db, model, param, true, verbose);

  /* Processing */

  for (int isimu = 0; isimu < getNbSimu(); isimu++)
  {

    /* Initiate the random normal values */

    _defineRandom();

    /* Apply the symmetry */

    _defineSymmetry();

    /* Perform the simulation */

    _final(db, iptr + isimu);
  }

  return 0;
}

/****************************************************************************/
/*!
 **  Dimension the ST_FFT structure
 **
 ** \param[in]  db      Db structure
 ** \param[in]  model   Model structure
 ** \param[in]  param   SimuFFTParam structure
 ** \param[in]  verbose Verbosity flag
 **
 *****************************************************************************/
void SimuFFT::_alloc(DbGrid *db,
                            Model *model,
                            const SimuFFTParam& param,
                            bool verbose)
{
  _nx.resize(3,0);
  _shift.resize(3,0);
  _dims.resize(3,0);
  _dim2.resize(3,0);
  _nxyz = 1;
  for (int i = 0; i < 3; i++)
  {
    if (i < _ndim)
    {
      _nx[i] = db->getNX(i);
      _nxyz *= _nx[i];
    }
  }

  /* Dilate the grid */

  _gridDilate(db, model, param, verbose);

  /* Determine the grid extension */

  for (int i = 0; i < 3; i++)
  {
    if (i < _ndim)
    {
      int nval = _getOptimalEvenNumber(_shift[i] + db->getNX(i));
      _dims[i] = nval;
      _dim2[i] = nval / 2;
    }
    else
    {
      _dims[i] = 1;
      _dim2[i] = 0;
    }
  }

  int total = 1;
  for (int idim = 0; idim < _ndim; idim++)
    total *= _dims[idim];
  _sizes_alloc = total;

  if (verbose)
  {
    message("Grid parameters after Optimal Dilation :\n");
    if (_ndim >= 1) message("- Number of Nodes along X = %d\n", _dims[0]);
    if (_ndim >= 2) message("- Number of Nodes along Y = %d\n", _dims[1]);
    if (_ndim >= 3) message("- Number of Nodes along Z = %d\n", _dims[2]);
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
int SimuFFT::_getOptimalEvenNumber(int number, int largeFactor)
{
  int local = number;
  if ((local % 2) == 1) local++;

  bool answer = true;
  while (answer)
  {
    VectorInt factors = _getFactors(local);
    int nfact = (int) factors.size();
    answer = false;
    for (int i = 0; i < nfact; i++)
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
VectorInt SimuFFT::_getFactors(int number)
{
  VectorInt factors;
  int local = number;
  int nfact = 0;

  /* Decomposition in multiples of 2 */

  int j = 2;
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
  }
  while (j <= local);

  if (nfact <= 0)
    factors.push_back(1);

  return factors;
}

/****************************************************************************/
/*!
 **  Calculates the grid extension in a given grid direction
 **
 ** \param[in]   db      Db structure
 ** \param[in]   model   Model structure
 ** \param[in]   param   SimuFFTParam structure
 ** \param[in]   verbose Verbosity flag
 **
 *****************************************************************************/
void SimuFFT::_gridDilate(const DbGrid *db,
                          Model *model,
                          const SimuFFTParam& param,
                          bool verbose)
{
  VectorInt indg(3);
  VectorDouble xyz0(3);
  VectorVectorDouble xyz(3);
  double percent = param.getPercent();

  /* Origin of the grid */

  for (int i = 0; i < 3; i++) indg[i] = 0;
  db->rankToCoordinate(db->indiceToRank(indg), xyz0);

  /* Location of the elementary end point */

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++) indg[j] = 0;
    indg[i] = 1;
    xyz[i].resize(3);
    db->rankToCoordinate(db->indiceToRank(indg), xyz[i]);
  }

  /* Coordinates of the grid vector in the rotated space */

  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
    {
      xyz[j][i] -= xyz0[i];
      if (j >= _ndim) xyz[j][i] = 0.;
    }

  /* Evaluate the count of elementary grid mesh (in each direction) */
  /* for the covariance to become negligeable (<percent)*/

  int ndx = 1;
  int ndy = 1;
  int ndz = 1;
  bool not_ok = true;

  while (not_ok)
  {
    not_ok = false;

    /* Extension along X */

    if (_ndim >= 1)
    {
      bool not_ok_dir = true;
      while (not_ok_dir)
      {
        bool correct = true;
        for (int idy = 0; idy < ndy && correct; idy++)
          for (int idz = 0; idz < ndz && correct; idz++)
            correct = _checkCorrect(model, xyz, ndx, idy, idz, percent);

        if (correct)
          not_ok_dir = false;
        else
        {
          ndx++;
          not_ok = true;
        }
      }
    }

    /* Extension along Y */

    if (_ndim >= 2)
    {
      bool not_ok_dir = true;
      while (not_ok_dir)
      {
        bool correct = true;
        for (int idx = 0; idx < ndx && correct; idx++)
          for (int idz = 0; idz < ndz && correct; idz++)
            correct = _checkCorrect(model, xyz, idx, ndy, idz, percent);

        if (correct)
          not_ok_dir = false;
        else
        {
          ndy++;
          not_ok = true;
        }
      }
    }

    /* Extension along Z */

    if (_ndim >= 3)
    {
      bool not_ok_dir = true;
      while (not_ok_dir)
      {
        bool correct = true;
        for (int idx = 0; idx < ndx && correct; idx++)
          for (int idy = 0; idy < ndy && correct; idy++)
            correct = _checkCorrect(model, xyz, idx, idy, ndz, percent);

        if (correct)
          not_ok_dir = false;
        else
        {
          ndz++;
          not_ok = true;
        }
      }
    }
  }

  /* Ultimate corrections */

  if (_ndim < 3) ndz = 0;
  if (_ndim < 2) ndy = 0;
  if (_ndim < 1) ndx = 0;

  /* Storing arguments */

  _shift[0] = ndx;
  _shift[1] = ndy;
  _shift[2] = ndz;

  /* Optional printout */

  if (verbose)
  {
    message("Grid Dilation parameters :\n");
    if (_ndim >= 1) message("- Number of Nodes along X = %d\n", ndx);
    if (_ndim >= 2) message("- Number of Nodes along Y = %d\n", ndy);
    if (_ndim >= 3) message("- Number of Nodes along Z = %d\n", ndz);
  }

  return;
}
/****************************************************************************/
/*!
 **  Checks if the covariance is below threshold for tested distance
 **
 ** \return True if the grid node is below threshold; 0 otherwise
 **
 ** \param[in]   model   Model structure
 ** \param[in]   xyz     Grid increment
 ** \param[in]   ix      Grid index along X
 ** \param[in]   iy      Grid index along Y
 ** \param[in]   iz      Grid index along Z
 ** \param[in]   percent Percentage of the model variance below which the
 **                      covariance is considered as small enough for dilation
 **
 *****************************************************************************/
bool SimuFFT::_checkCorrect(Model *model,
                            const VectorVectorDouble& xyz,
                            int ix,
                            int iy,
                            int iz,
                            double percent)
{
  double hh, value, refval;
  VectorDouble d(3);

  /* Calculate the reference C(0) value */

  for (int i = 0; i < 3; i++) d[i] = 0.;
  hh = ut_vector_norm(d);
  (void) model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0, ECalcMember::LHS, 1, d,
                        &hh, &refval);

  /* Calculate the distance */

  for (int i = 0; i < 3; i++)
    d[i] = ix * xyz[i][0] + iy * xyz[i][1] + iz * xyz[i][2];
  hh = ut_vector_norm(d);

  /* Evaluate the covariance value */

  (void) model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0, ECalcMember::LHS, 1, d,
                        &hh, &value);

  if (value / refval > percent / 100) return false;
  return true;
}

/****************************************************************************/
/*!
 **  Calculate the norme of a vector
 **
 ** \return  Norm of the vector
 **
 ** \param[in]  incr  Increment vector
 **
 *****************************************************************************/
double st_norm(VectorDouble &incr)

{
  ut_vector_norm(incr);
  int i;
  double value;

  value = 0.;
  for (i = 0; i < 3; i++)
    value += incr[i] * incr[i];
  value = sqrt(value);
  return (value);
}

/****************************************************************************/
/*!
 **  Prepares the simulation on a grid using Discrete FFT
 **
 ** \param[in]  db              Db structure
 ** \param[in]  model           Model structure
 ** \param[in]  param           SimuFFTParam structure
 ** \param[in]  flag_amplitude  1 to convert into amplitude
 ** \param[in]  verbose         Verbosity flag
 ** \param[in]  eps             Tolerance
 **
 *****************************************************************************/
void SimuFFT::_prepar(DbGrid *db,
                      Model *model,
                      const SimuFFTParam& param,
                      bool flag_amplitude,
                      bool verbose,
                      double eps)
{
  double value;
  VectorDouble del(3);
  VectorDouble xyz0(3);
  VectorDouble xyz(3);
  VectorDouble delta(3);
  VectorInt indg(3);
  VectorInt jnd(3);
  VectorVectorDouble xyz1(3);

  /* Initializations */

  int kbmax = (param.isFlagAliasing()) ? 1 : 0;
  VectorDouble cplx(_sizes_alloc);
  VectorDouble cply(_sizes_alloc);

  /* Local core allocation */

  double hnorm = 1.;
  for (int i = 0; i < 3; i++)
  {
    indg[i] = 0;
    delta[i] = del[i] = 0.;
    xyz[i] = xyz0[i] = 0.;
    xyz1[i].resize(3);
    for (int j = 0; j < 3; j++)
      xyz1[i][j] = 0.;
  }
  db->rankToCoordinate(db->indiceToRank(indg), xyz0);

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++) indg[j] = 0;
    indg[i] = 1;
    xyz1[i].resize(3);
    if (i < _ndim)
    {
      db->rankToCoordinate(db->indiceToRank(indg), xyz1[i]);
      for (int j = 0; j < 3; j++)
        xyz1[i][j] -= xyz0[j];
      delta[i] = db->getDX(i) * db->getNX(i);
    }
    else
    {
      for (int j = 0; j < 3; j++)
        xyz1[i][j] = 0.;
      delta[i] = 0.;
    }
  }

  /* Loop for anti-aliasing */

  for (int kbound = 0; kbound <= kbmax; kbound++)
  {

    /* Local initializations */

    int kb1 = (_ndim >= 1) ? kbound : 0;
    int kb2 = (_ndim >= 2) ? kbound : 0;
    int kb3 = (_ndim >= 3) ? kbound : 0;
    for (int i = 0; i < _sizes_alloc; i++)
      cplx[i] = cply[i] = 0.;

    /* Calculate the normation scale */

    double scale = 0.;
    for (int k1 = -kb1; k1 <= kb1; k1++)
      for (int k2 = -kb2; k2 <= kb2; k2++)
        for (int k3 = -kb3; k3 <= kb3; k3++)
        {
          del[0] = k1 * delta[0];
          del[1] = k2 * delta[1];
          del[2] = k3 * delta[2];
          hnorm = ut_vector_norm(del);
          (void) model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0,
                                ECalcMember::LHS, 1, del, &hnorm, &value);
          scale += value;
        }
    for (int i = 0; i < 3; i++)
      del[i] = 0.;
    hnorm = ut_vector_norm(del);
    (void) model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0, ECalcMember::LHS, 1,
                          del, &hnorm, &value);
    double coeff = value / scale;

    int ecr = 0;
    for (int iz = 0; iz < _dims[2]; iz++)
      for (int iy = 0; iy < _dims[1]; iy++)
        for (int ix = 0; ix < _dims[0]; ix++, ecr++)
        {
          jnd[0] = (ix <= _dim2[0]) ? ix : ix - _dims[0];
          jnd[1] = (iy <= _dim2[1]) ? iy : iy - _dims[1];
          jnd[2] = (iz <= _dim2[2]) ? iz : iz - _dims[2];
          for (int i = 0; i < 3; i++)
          {
            double proj = 0.;
            for (int j = 0; j < 3; j++)
              proj += jnd[j] * xyz1[j][i];
            xyz[i] = (i < _ndim) ? proj : 0.;
          }
          for (int k1 = -kb1; k1 <= kb1; k1++)
            for (int k2 = -kb2; k2 <= kb2; k2++)
              for (int k3 = -kb3; k3 <= kb3; k3++)
              {
                del[0] = xyz[0] + k1 * delta[0];
                del[1] = xyz[1] + k2 * delta[1];
                del[2] = xyz[2] + k3 * delta[2];
                hnorm = st_norm(del);
                (void) model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0,
                                      ECalcMember::LHS, 1, del, &hnorm, &value);
                cplx[ecr] += coeff * value;
              }
        }

    /* Perform the Fast Fourier Transform */

    (void) fftn(_ndim, _dims.data(), cplx.data(), cply.data(), -1, 1.);

    /* Looking for negative terms */

    double total_plus = 0.;
    double total_moins = 0.;
    for (int i = 0; i < _sizes_alloc; i++)
    {
      cplx[i] /= (double) _sizes_alloc;
      if (cplx[i] < 0)
        total_moins -= cplx[i];
      else
        total_plus += cplx[i];
    }

    /* Correcting positive terms of the spectrum */

    double correc = (total_plus - total_moins) / total_plus;
    if (total_moins > 0)
    {
      for (int i = 0; i < _sizes_alloc; i++)
        if (cplx[i] < 0.)
          cplx[i] = 0.;
        else
          cplx[i] *= correc;
    }

    /* Converting into amplitude */

    if (flag_amplitude)
      for (int i = 0; i < _sizes_alloc; i++)
        cplx[i] = sqrt(cplx[i] / 2.);

    // Returned arguments

    _cmat = cplx;
    _rnd  = cply;

    /* Printout statistics */

    if (verbose)
    {
      message("Statistics on the Discrete Periodic Covariance\n");
      if (param.isFlagAliasing())
        message("- Anti-aliasing switched ON (Iteration=%d)\n", kbound + 1);
      else
        message("- Anti-aliasing switched OFF\n");
      message("- Total of positive frequencies = %lf\n", total_plus);
      message("- Total of negative frequencies = %lf\n", total_moins);
      message("\n");
    }
    if (ABS(correc - 1.) < eps) break;
  }

  return;
}

/****************************************************************************/
/*!
 **  Initiate a vector of random normal values
 **
 *****************************************************************************/
void SimuFFT::_defineRandom()

{
  for (int i = 0; i < _sizes_alloc; i++)
    _u[i] = _cmat[i] * law_gaussian();
  for (int i = 0; i < _sizes_alloc; i++)
    _v[i] = _cmat[i] * law_gaussian();

  switch (_ndim)
  {
    case 1:
      for (int ix = 0; ix < _dims[0]; ix += _dim2[0])
        _setVariance(ix, 0, 0);
      break;

    case 2:
      for (int iy = 0; iy < _dims[1]; iy += _dim2[1])
        for (int ix = 0; ix < _dims[0]; ix += _dim2[0])
          _setVariance(ix, iy, 0);
      break;

    case 3:
      for (int iz = 0; iz < _dims[2]; iz += _dim2[2])
        for (int iy = 0; iy < _dims[1]; iy += _dim2[1])
          for (int ix = 0; ix < _dims[0]; ix += _dim2[0])
            _setVariance(ix, iy, iz);
      break;
  }
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
void SimuFFT::_setVariance(int ix, int iy, int iz)
{
  int ind = IND(ix, iy, iz);
  _u[ind] *= sqrt(2.0);
  _v[ind] = 0.;
}

/****************************************************************************/
/*!
 **  Operate the symmetry
 **
 *****************************************************************************/
void SimuFFT::_defineSymmetry(void)
{

  /* Dispatch according to the space dimension */

  switch (_ndim)
  {
    case 1:
      _defineSym1();
      break;

    case 2:
      _defineSym2(0);
      break;

    case 3:
      _defineSym3();
      break;
  }
  return;
}

/****************************************************************************/
/*!
 **  Operate the symmetry for a 1-D space
 **
 *****************************************************************************/
void SimuFFT::_defineSym1()

{
  // A(1) and A(N1/2+1) are real
  for (int ix = 0; ix < _dims[0]; ix += _dim2[0])
    _setZero(ix, 0, 0);

  // A(j) = A*(N1-j+2) for j in [2,N1/2]
  for (int ix = 1; ix < _dim2[0]; ix++)
  {
    int jx = _dims[0] - ix;
    _setConjugate(ix, 0, 0, jx, 0, 0);
  }

  return;
}

/****************************************************************************/
/*!
 **  Operate the symmetry for a 2-D space
 **
 ** \param[in]  iz0   fixed third index
 **
 *****************************************************************************/
void SimuFFT::_defineSym2(int iz0)
{
  // A(1,1), A(N1/2,1), A(1,N2/2) and A(N1/2,N2/2) real
  for (int iy = 0; iy < _dims[1]; iy += _dim2[1])
    for (int ix = 0; ix < _dims[0]; ix += _dim2[0])
      _setZero(ix, iy, iz0);

  // A(1,k)      = A*(1,N2-k+2)      for k in [2,N2/2]
  // A(N2/2+1,k) = A*(N2/2+1,N2-k+2) for k in [2,N2/2]
  for (int iy = 1; iy < _dim2[1]; iy++)
    for (int ix = 0; ix < _dims[0]; ix += _dim2[0])
    {
      int jy = _dims[1] - iy;
      _setConjugate(ix, iy, iz0, ix, jy, iz0);
    }

  // A(j,1)      = A*(N1-j+2,1)      for j in [2,N1/2]
  // A(j,N2/2+1) = A*(N1-j+2,N2/2+1) for j in [2,N1/2]
  for (int iy = 0; iy < _dims[1]; iy += _dim2[1])
    for (int ix = 1; ix < _dim2[0]; ix++)
    {
      int jx = _dims[0] - ix;
      _setConjugate(ix, iy, iz0, jx, iy, iz0);
    }

  // A(j,k) = A*(N1-j+2,N2-k+2) for j in [2,N1/2] and k in [2,N2/2]
  for (int iy = 1; iy < _dim2[1]; iy++)
    for (int ix = 1; ix < _dim2[0]; ix++)
    {
      int jx = _dims[0] - ix;
      int jy = _dims[1] - iy;
      _setConjugate(ix, iy, iz0, jx, jy, iz0);
    }

  // A(j,N2-k+2) = A*(N1-j+2,k) for j in [2,N1/2] and k in [2,N2/2]
  for (int iy = 1; iy < _dim2[1]; iy++)
    for (int ix = 1; ix < _dim2[0]; ix++)
    {
      int jx = _dims[0] - ix;
      int jy = _dims[1] - iy;
      _setConjugate(ix, jy, iz0, jx, iy, iz0);
    }

  return;
}

/****************************************************************************/
/*!
 **  Operate the symmetry for a 3-D space
 **
 *****************************************************************************/
void SimuFFT::_defineSym3()

{
  // For l=1 or N3/2+1, use the 2-D symmetry
  for (int iz = 0; iz < _dims[2]; iz += _dim2[2])
    _defineSym2(iz);

  // For the other planes:

  // A(1,1,l)           = A*(1,1,N3-l+2)           for l in [2,N3/2]
  // A(N1/2+1,1,l)      = A*(N1/2+1,1,N3-l+2)      for l in [2,N3/2]
  // A(1,N2/2+1,l)      = A*(1,N2/2+1,N3-l+2)      for l in [2,N3/2]
  // A(N1/2+1,N2/2+1,l) = A*(N1/2+1,N2/2+1,N3-l+2) for l in [2,N3/2]
  for (int iz = 1; iz < _dim2[2]; iz++)
    for (int iy = 0; iy < _dims[1]; iy += _dim2[1])
      for (int ix = 0; ix < _dims[0]; ix += _dim2[0])
      {
        int jz = _dims[2] - iz;
        _setConjugate(ix, iy, iz, ix, iy, jz);
      }

  // A(1,k,l)      = A*(1,N2-k+2,N3-l+2)      for k in [2,N2/2] and l in [2,N3/2]
  // A(N1/2+1,k,l) = A*(N1/2+1,N2-k+2,N3-l+2) for k in [2,N2/2] and l in [2,N3/2]
  for (int iz = 1; iz < _dim2[2]; iz++)
    for (int iy = 1; iy < _dim2[1]; iy++)
      for (int ix = 0; ix < _dims[0]; ix += _dim2[0])
      {
        int jy = _dims[1] - iy;
        int jz = _dims[2] - iz;
        _setConjugate(ix, iy, iz, ix, jy, jz);
      }

  // A(1,N2-k+2,l)      = A*(1,k,N3/-l+2)      for k in [2,N2/2] and l in [2,N3/2]
  // A(N1/2+1,N2-k+2,l) = A*(N1/2+1,k,N3/-l+2) for k in [2,N2/2] and l in [2,N3/2]
  for (int iz = 1; iz < _dim2[2]; iz++)
    for (int iy = 1; iy < _dim2[1]; iy++)
      for (int ix = 0; ix < _dims[0]; ix += _dim2[0])
      {
        int jy = _dims[1] - iy;
        int jz = _dims[2] - iz;
        _setConjugate(ix, jy, iz, ix, iy, jz);
      }

  // A(j,1,l)      = A*(N1-j+2,1,N3-l+2)      for j in [2,N1/2] and l in [2,N3/2]
  // A(j,N2/2+1,l) = A*(N1-j+2,N2/2+1,N3-l+2) for j in [2,N1/2] and l in [2,N3/2]
  for (int iz = 1; iz < _dim2[2]; iz++)
    for (int iy = 0; iy < _dims[1]; iy += _dim2[1])
      for (int ix = 1; ix < _dim2[0]; ix++)
      {
        int jx = _dims[0] - ix;
        int jz = _dims[2] - iz;
        _setConjugate(ix, iy, iz, jx, iy, jz);
      }

  // A(N1-j+2,1,l)      = A*(j,1,N3-l+2)      for j in [2,N1/2] and l in [2,N3/2]
  // A(N1-j+2,N2/2+1,l) = A*(j,N2/2+1,N3-l+2) for j in [2,N1/2] and l in [2,N3/2]
  for (int iz = 1; iz < _dim2[2]; iz++)
    for (int iy = 0; iy < _dims[1]; iy += _dim2[1])
      for (int ix = 1; ix < _dim2[0]; ix++)
      {
        int jx = _dims[0] - ix;
        int jz = _dims[2] - iz;
        _setConjugate(jx, iy, iz, ix, iy, jz);
      }

  // A(j,k,l) = A*(N1-j+2,N2-k+2,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
  for (int iz = 1; iz < _dim2[2]; iz++)
    for (int iy = 1; iy < _dim2[1]; iy++)
      for (int ix = 1; ix < _dim2[0]; ix++)
      {
        int jx = _dims[0] - ix;
        int jy = _dims[1] - iy;
        int jz = _dims[2] - iz;
        _setConjugate(ix, iy, iz, jx, jy, jz);
      }

  // A(N1-j+2,N2-k+2,l) = A*(j,k,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
  for (int iz = 1; iz < _dim2[2]; iz++)
    for (int iy = 1; iy < _dim2[1]; iy++)
      for (int ix = 1; ix < _dim2[0]; ix++)
      {
        int jx = _dims[0] - ix;
        int jy = _dims[1] - iy;
        int jz = _dims[2] - iz;
        _setConjugate(jx, jy, iz, ix, iy, jz);
      }

  // A(j,N2-k+2,l) = A*(N2-j+2,k,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
  for (int iz = 1; iz < _dim2[2]; iz++)
    for (int iy = 1; iy < _dim2[1]; iy++)
      for (int ix = 1; ix < _dim2[0]; ix++)
      {
        int jx = _dims[0] - ix;
        int jy = _dims[1] - iy;
        int jz = _dims[2] - iz;
        _setConjugate(ix, jy, iz, jx, iy, jz);
      }

  // A(N1-j+2,k,l) = A*(j,k,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
  for (int iz = 1; iz < _dim2[2]; iz++)
    for (int iy = 1; iy < _dim2[1]; iy++)
      for (int ix = 1; ix < _dim2[0]; ix++)
      {
        int jx = _dims[0] - ix;
        int jy = _dims[1] - iy;
        int jz = _dims[2] - iz;
        _setConjugate(jx, iy, iz, ix, jy, jz);
      }

  return;
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
void SimuFFT::_setZero(int ix, int iy, int iz)
{
  int ind = IND(ix, iy, iz);
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
void SimuFFT::_setConjugate(int ix, int iy, int iz, int jx, int jy, int jz)
{
  int ind1 = IND(ix, iy, iz);
  int ind2 = IND(jx, jy, jz);
  _u[ind2] = _u[ind1];
  _v[ind2] = -_v[ind1];
}

/****************************************************************************/
/*!
 **  Perform a non-conditional simulation on the grid
 **
 ** \param[in]  db    Db structure
 ** \param[in]  iad   address for writing the simulation
 **
 *****************************************************************************/
void SimuFFT::_final(DbGrid *db, int iad)
{
  (void) fftn(_ndim, _dims.data(), _u.data(), _v.data(), 1, 1.);
  int nx = MAX(_nx[0], 1);
  int ny = MAX(_nx[1], 1);
  int nz = MAX(_nx[2], 1);

  /* Retrieving the simulation */

  int ecr = 0;
  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++, ecr++)
      {
        int jx = ix + _shift[0];
        int jy = iy + _shift[1];
        int jz = iz + _shift[2];
        db->updArray(ecr, iad, 4, U(jx, jy, jz));
      }

  return;
}

/****************************************************************************/
/*!
 **  Calculate the change of support coefficients by FFT method
 **  in the lognormal case on a grid
 **
 ** \return  r^2 coefficients for the different logarithmic variances
 **
 ** \param[in]  db      Db structure
 ** \param[in]  model   Model structure
 ** \param[in]  param   SimuFFTParam structure
 ** \param[in]  sigma   Array of logarithmic variances
 ** \param[in]  verbose Verbosity flag
 **
 *****************************************************************************/
VectorDouble SimuFFT::getChangeSupport(DbGrid *db,
                                       Model *model,
                                       const SimuFFTParam& param,
                                       const VectorDouble& sigma,
                                       bool verbose)
{
  // Check validity of the context

  if (! _isValid(db, model)) return VectorDouble();

  _alloc(db, model, param, verbose);

  /* Preparation of the FFT environment */

  _prepar(db, model, param, false, verbose);

  /* Calculate the correlation matrix (possibly rescaled) */

  (void) fftn(_ndim, _dims.data(), _cmat.data(), _rnd.data(), 1, 1.);

  /* Loop on the different lognormal variances */

  int nval = (int) sigma.size();
  VectorDouble r2val;
  if (nval > 0)
  {
    r2val.resize(nval,TEST);
    for (int ival = 0; ival < nval; ival++)
      r2val[ival] = _support(sigma[ival]);
  }
  else
  {
    r2val.resize(1,TEST);
    r2val[0] = _support(TEST);
  }
  return r2val;
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
double SimuFFT::_support(double sigma)
{
  double value = 0.;
  if (sigma == 0.) return (TEST);

  switch (_ndim)
  {
    case 1:
      value = _support1(sigma);
      break;

    case 2:
      value = _support2(sigma);
      break;

    case 3:
      value = _support3(sigma);
      break;
  }

  /* Calculate the scale */

  double scale = 1.;
  for (int idim = 0; idim < _ndim; idim++)
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
double SimuFFT::_support1(double sigma)
{
  double value = 0.;
  for (int ix = -_nx[0]; ix <= _nx[0]; ix++)
  {
    int iix = (ix < 0) ? _dims[0] + ix : ix;
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
double SimuFFT::_support2(double sigma)
{
  double value = 0.;
  for (int ix = -_nx[0]; ix <= _nx[0]; ix++)
    for (int iy = -_nx[1]; iy <= _nx[1]; iy++)
    {
      int iix = (ix < 0) ? _dims[0] + ix : ix;
      int iiy = (iy < 0) ? _dims[1] + iy : iy;
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
double SimuFFT::_support3(double sigma)
{
  double value = 0.;
  for (int ix = -_nx[0]; ix <= _nx[0]; ix++)
    for (int iy = -_nx[1]; iy <= _nx[1]; iy++)
      for (int iz = -_nx[2]; iz <= _nx[2]; iz++)
      {
        int iix = (ix < 0) ? _dims[0] + ix : ix;
        int iiy = (iy < 0) ? _dims[1] + iy : iy;
        int iiz = (iz < 0) ? _dims[2] + iz : iz;
        double rho = _rhoSigma(sigma, iix, iiy, iiz);
        value += ((_nx[0] - ABS(ix)) * (_nx[1] - ABS(iy)) * (_nx[2] - ABS(iz)) * rho);
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
double SimuFFT::_rhoSigma(double sigma, int ix, int iy, int iz)
{
  double rho = _cmat[IND(ix, iy, iz)];
  if (!FFFF(sigma)) rho = exp(sigma * sigma * rho);
  return (rho);
}

/****************************************************************************/
/*!
 **  Checks the environment for simulations
 **
 ** \return  Error return code
 **
 ** \param[in]  db     grid Db structure
 ** \param[in]  model  Model structure
 **
 *****************************************************************************/
bool SimuFFT::_isValid(Db *db, Model *model)
{
  _ndim = db->getNDim();

  /**************************************************************/
  /* Check if the Space dimension is compatible with the method */
  /**************************************************************/

  if (_ndim < 1 || _ndim > 3)
  {
    messerr("The FFT Method is not a relevant simulation model");
    messerr("for this Space Dimension (%d)", _ndim);
    return false;
  }

  /**********************/
  /* Checking the model */
  /**********************/

  if (model != nullptr)
  {
    if (model->getVariableNumber() != 1)
    {
      messerr("The FFT method is restricted to the monovariate case (%d)",
              model->getVariableNumber());
      return false;
    }
    if (model->getCovaNumber() <= 0)
    {
      messerr("The number of covariance must be positive");
      return false;
    }
    if (model->getDimensionNumber() <= 0)
    {
      messerr("The Space Dimension must be positive = %d",
              model->getDimensionNumber());
      return false;
    }
    if (model->getDimensionNumber() != _ndim)
    {
      messerr("The Space Dimension of the Db structure (%d)", _ndim);
      messerr("Does not correspond to the Space Dimension of the model (%d)",
              model->getDimensionNumber());
      return false;
    }
  }

  /********************/
  /* Calculate the Db */
  /********************/

  VectorDouble db_mini(_ndim);
  VectorDouble db_maxi(_ndim);
  db_extension(db, db_mini, db_maxi);

  if (model != nullptr)
    model->setField(ut_vector_extension_diagonal(db_mini, db_maxi));

  return true;
}

bool SimuFFT::_run()
{
  return true;
}
