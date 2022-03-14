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
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"

#include <math.h>

/*! \cond */
#define LARGE_FACTOR 11
#define MAXFACTOR 100
#define IND(ix,iy,iz) ((iz) + simu->dims[2] * ((iy) + simu->dims[1] * (ix)))
#define U(ix,iy,iz)   (simu->u[IND(ix,iy,iz)])
#define V(ix,iy,iz)   (simu->v[IND(ix,iy,iz)])
#define EPS 1.e-5
#define DEBUG 0
/*! \endcond */

typedef struct
{
  int ndim;
  int nxyz;
  int nx[3];
  int shift[3];
  int dims[3];
  int dim2[3];
  int sizes_alloc;
  double *cmat;
  double *rnd;
  double *u;
  double *v;
} ST_FFT;

static int FACTORS[MAXFACTOR];
static int FLAG_ALIASING = 0;

/****************************************************************************/
/*!
 **  Calculate the norme of a vector
 **
 ** \return  Norm of the vector
 **
 ** \param[in]  incr  Increment vector
 **
 *****************************************************************************/
static double st_norm(VectorDouble &incr)

{
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
 **  Get the factor decomposition of a number
 **
 ** \return  Count of active factors
 **
 ** \param[in]  number  number to be decomposed
 **
 *****************************************************************************/
static int st_get_factors(int number)

{
  int j, local, nfact;

  /* Initializations */

  local = number;
  nfact = 0;

  /* Decomposition in multiples of 2 */

  j = 2;
  while ((local % j) == 0)
  {
    FACTORS[nfact++] = j;
    if (nfact >= MAXFACTOR) messageAbort("st_get_factors");
    local /= j;
  }

  /* Decomposition in higher level multiples */

  j = 3;
  do
  {
    while ((local % j) == 0)
    {
      FACTORS[nfact++] = j;
      if (nfact >= MAXFACTOR) messageAbort("st_get_factors");
      local /= j;
    }
    j += 2;
  }
  while (j <= local);

  if (nfact <= 0) FACTORS[nfact++] = 1;

  return (nfact);
}

/****************************************************************************/
/*!
 **  Returns the closest value, larger than the argument, which is
 **  factorized as the product of low factors
 **
 ** \return  Returned number
 **
 ** \param[in]  number input number
 **
 *****************************************************************************/
static int st_get_optimal_even_number(int number)

{
  int i, nfact, local, answer;

  local = number;
  if ((local % 2) == 1) local++;

  answer = 1;
  while (answer)
  {
    nfact = st_get_factors(local);
    for (i = answer = 0; i < nfact; i++)
      if (FACTORS[i] > LARGE_FACTOR) answer = 1;
    if (answer) local += 2;
  }
  return (local);
}

/****************************************************************************/
/*!
 **  Initialize the ST_FFT structure
 **
 ** \param[in]  db    Db description
 **
 ** \param[out] simu  Initialized structure
 **
 *****************************************************************************/
static void st_simfft_init(DbGrid *db, ST_FFT *simu)
{
  int i;

  /* Initializations */

  simu->sizes_alloc = 0;
  for (i = 0; i < 3; i++)
  {
    simu->shift[i] = 0;
    simu->dims[i] = 0;
    simu->dim2[i] = 0;
  }
  simu->cmat = nullptr;
  simu->rnd = nullptr;
  simu->u = nullptr;
  simu->v = nullptr;

  /* Definition according to the grid */

  simu->ndim = db->getNDim();
  simu->nxyz = 1;
  for (i = 0; i < 3; i++)
  {
    simu->nx[i] = db->getNX(i);
    simu->nxyz *= simu->nx[i];
  }

  return;
}

/****************************************************************************/
/*!
 **  Prepares the simulation on a grid using Discrete FFT
 **
 ** \param[in]  db              Db structure
 ** \param[in]  model           Model structure
 ** \param[in]  simu            ST_FFT structure
 ** \param[in]  flag_amplitude  1 to convert into amplitude
 **
 *****************************************************************************/
static void st_simfft_prepar(DbGrid *db,
                             Model *model,
                             ST_FFT *simu,
                             int flag_amplitude)
{
  double *cplx, *cply, proj, xyz0[3], xyz1[3][3], xyz[3], delta[3];
  double total_plus, total_moins, correc, scale, coeff, hnorm, value;
  int i, j, k1, k2, k3, ix, iy, iz, jnd[3], ecr, indg[3];
  int kbound, kbmax, kb1, kb2, kb3, ndim;
  VectorDouble del;

  /* Initializations */

  kbmax = (FLAG_ALIASING) ? 1 :
                            0;
  ndim = simu->ndim;
  cplx = simu->cmat;
  cply = simu->rnd;
  del.resize(3);

  /* Local core allocation */

  hnorm = 1.;
  for (i = 0; i < 3; i++)
  {
    indg[i] = jnd[i] = 0;
    delta[i] = del[i] = 0.;
    xyz[i] = xyz0[i] = 0.;
    for (j = 0; j < 3; j++)
      xyz1[i][j] = 0.;
  }
  grid_to_point(db, indg, nullptr, xyz0);
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
      indg[j] = 0;
    indg[i] = 1;
    grid_to_point(db, indg, nullptr, xyz1[i]);
    for (j = 0; j < 3; j++)
      xyz1[i][j] -= xyz0[j];
    delta[i] = db->getDX(i) * db->getNX(i);
    if (i < ndim) continue;
    delta[i] = 0.;
    for (j = 0; j < 3; j++)
      xyz1[i][j] = 0.;
  }

  /* Loop for anti-aliasing */

  for (kbound = 0; kbound <= kbmax; kbound++)
  {

    /* Local initializations */

    kb1 = (ndim >= 1) ? kbound :
                        0;
    kb2 = (ndim >= 2) ? kbound :
                        0;
    kb3 = (ndim >= 3) ? kbound :
                        0;
    for (i = 0; i < simu->sizes_alloc; i++)
      cplx[i] = cply[i] = 0.;

    /* Calculate the normation scale */

    scale = 0.;
    for (k1 = -kb1; k1 <= kb1; k1++)
      for (k2 = -kb2; k2 <= kb2; k2++)
        for (k3 = -kb3; k3 <= kb3; k3++)
        {
          del[0] = k1 * delta[0];
          del[1] = k2 * delta[1];
          del[2] = k3 * delta[2];
          hnorm = st_norm(del);
          (void) model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0,
                                ECalcMember::LHS, 1, del, &hnorm, &value);
          scale += value;
        }
    for (i = 0; i < 3; i++)
      del[i] = 0.;
    hnorm = st_norm(del);
    (void) model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0, ECalcMember::LHS, 1,
                          del, &hnorm, &value);
    coeff = value / scale;

    ecr = 0;
    for (iz = 0; iz < simu->dims[2]; iz++)
      for (iy = 0; iy < simu->dims[1]; iy++)
        for (ix = 0; ix < simu->dims[0]; ix++, ecr++)
        {
          jnd[0] = (ix <= simu->dim2[0]) ? ix :
                                           ix - simu->dims[0];
          jnd[1] = (iy <= simu->dim2[1]) ? iy :
                                           iy - simu->dims[1];
          jnd[2] = (iz <= simu->dim2[2]) ? iz :
                                           iz - simu->dims[2];
          for (i = 0; i < 3; i++)
          {
            proj = 0.;
            for (j = 0; j < 3; j++)
              proj += jnd[j] * xyz1[j][i];
            xyz[i] = (i < ndim) ? proj :
                                  0.;
          }
          for (k1 = -kb1; k1 <= kb1; k1++)
            for (k2 = -kb2; k2 <= kb2; k2++)
              for (k3 = -kb3; k3 <= kb3; k3++)
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

    if (DEBUG)
      print_matrix("Discretized Covariance", 0, 1, simu->dims[0], simu->dims[1],
                   NULL, cplx);
    (void) fftn(simu->ndim, simu->dims, cplx, cply, -1, 1.);

    /* Looking for negative terms */

    total_plus = total_moins = 0.;
    for (i = 0; i < simu->sizes_alloc; i++)
    {
      cplx[i] /= (double) simu->sizes_alloc;
      if (cplx[i] < 0)
        total_moins -= cplx[i];
      else
        total_plus += cplx[i];
    }

    /* Correcting positive terms of the spectrum */

    correc = (total_plus - total_moins) / total_plus;
    if (total_moins > 0)
    {
      for (i = 0; i < simu->sizes_alloc; i++)
        if (cplx[i] < 0.)
          cplx[i] = 0.;
        else
          cplx[i] *= correc;
    }

    /* Converting into amplitude */

    if (flag_amplitude) for (i = 0; i < simu->sizes_alloc; i++)
      cplx[i] = sqrt(cplx[i] / 2.);

    /* Printout statistics */

    if (OptDbg::query(EDbg::SIMULATE))
    {
      message("Statistics on the Discrete Periodic Covariance\n");
      if (FLAG_ALIASING)
        message("- Anti-aliasing switched ON (Iteration=%d)\n", kbound + 1);
      else
        message("- Anti-aliasing switched OFF\n");
      message("- Total of positive frequencies = %lf\n", total_plus);
      message("- Total of negative frequencies = %lf\n", total_moins);
      message("\n");
    }
    if (DEBUG)
      print_matrix("Discretized (square root of half-) Spectrum", 0, 1,
                   simu->dims[0], simu->dims[1], NULL, cplx);
    if (ABS(correc - 1.) < EPS) break;
  }

  return;
}

/****************************************************************************/
/*!
 **  Correct the variance of the spectrum for real U
 **
 ** \param[in]  simu  ST_FFT structure
 ** \param[in]  ix    Cell location along X
 ** \param[in]  iy    Cell location along Y
 ** \param[in]  iz    Cell location along Z
 **
 *****************************************************************************/
static void st_set_variance(ST_FFT *simu, int ix, int iy, int iz)
{
  int ind;

  ind = IND(ix, iy, iz);
  simu->u[ind] *= sqrt(2.0);
  simu->v[ind] = 0.;
  if (DEBUG) message("Element (%d,%d,%d) variance updated\n", ix, iy, iz);
}

/****************************************************************************/
/*!
 **  Set the imaginary part of a cell to zero
 **
 ** \param[in]  simu  ST_FFT structure
 ** \param[in]  ix    Cell location along X
 ** \param[in]  iy    Cell location along Y
 ** \param[in]  iz    Cell location along Z
 **
 *****************************************************************************/
static void st_set_zero(ST_FFT *simu, int ix, int iy, int iz)
{
  int ind;

  ind = IND(ix, iy, iz);
  simu->v[ind] = 0.;
  if (DEBUG) message("Element (%d,%d,%d) has no imaginary part\n", ix, iy, iz);
}

/****************************************************************************/
/*!
 **  Set the target cell as the conjugate of the input cell
 **
 ** \param[in]  simu  ST_FFT structure
 ** \param[in]  ix    Input cell location along X
 ** \param[in]  iy    Input cell location along Y
 ** \param[in]  iz    Input cell location along Z
 ** \param[in]  jx    Target cell location along X
 ** \param[in]  jy    Target cell location along Y
 ** \param[in]  jz    Target cell location along Z
 **
 *****************************************************************************/
static void st_set_conj(ST_FFT *simu,
                        int ix,
                        int iy,
                        int iz,
                        int jx,
                        int jy,
                        int jz)
{
  int ind1, ind2;

  ind1 = IND(ix, iy, iz);
  ind2 = IND(jx, jy, jz);
  simu->u[ind2] = simu->u[ind1];
  simu->v[ind2] = -simu->v[ind1];
  if (DEBUG)
    message("Copy the conjugate of (%d,%d,%d) into (%d,%d,%d)\n", ix, iy, iz,
            jx, jy, jz);
}

/****************************************************************************/
/*!
 **  Initiate a vector of random normal values
 **
 ** \param[in]  simu  ST_FFT structure
 **
 *****************************************************************************/
static void st_simfft_random(ST_FFT *simu)

{
  int i, ix, iy, iz;

  for (i = 0; i < simu->sizes_alloc; i++)
    simu->u[i] = simu->cmat[i] * law_gaussian();
  for (i = 0; i < simu->sizes_alloc; i++)
    simu->v[i] = simu->cmat[i] * law_gaussian();

  switch (simu->ndim)
  {
    case 1:
      for (ix = 0; ix < simu->dims[0]; ix += simu->dim2[0])
        st_set_variance(simu, ix, 0, 0);
      break;

    case 2:
      for (iy = 0; iy < simu->dims[1]; iy += simu->dim2[1])
        for (ix = 0; ix < simu->dims[0]; ix += simu->dim2[0])
          st_set_variance(simu, ix, iy, 0);
      break;

    case 3:
      for (iz = 0; iz < simu->dims[2]; iz += simu->dim2[2])
        for (iy = 0; iy < simu->dims[1]; iy += simu->dim2[1])
          for (ix = 0; ix < simu->dims[0]; ix += simu->dim2[0])
            st_set_variance(simu, ix, iy, iz);
      break;
  }
}

/****************************************************************************/
/*!
 **  Checks if the covariance is below threshold for tested distance
 **
 ** \param[in]   model   Model structure
 ** \param[in]   xyz     Grid increment
 ** \param[in]   ix      Grid index along X
 ** \param[in]   iy      Grid index along Y
 ** \param[in]   iz      Grid index along Z
 ** \param[in]   percent Percentage of the model variance below which the
 **                      covariance is considered as small enough for dilation
 **
 ** \param[out]  correct  1 if the grid node is below threshold; 0 otherwise
 **
 *****************************************************************************/
static void st_check_correct(Model *model,
                             double xyz[3][3],
                             int ix,
                             int iy,
                             int iz,
                             double percent,
                             int *correct)
{
  double hh, refval, value;
  int i;
  VectorDouble d;

  /* Calculate the reference C(0) value */

  d.resize(3);
  for (i = 0; i < 3; i++)
    d[i] = 0.;
  hh = st_norm(d);
  (void) model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0, ECalcMember::LHS, 1, d,
                        &hh, &refval);

  /* Calculate the distance */

  for (i = 0; i < 3; i++)
    d[i] = ix * xyz[i][0] + iy * xyz[i][1] + iz * xyz[i][2];
  hh = st_norm(d);

  /* Evaluate the covariance value */

  (void) model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0, ECalcMember::LHS, 1, d,
                        &hh, &value);

  if (value / refval > percent / 100) (*correct) = 0;

  return;
}

/****************************************************************************/
/*!
 **  Returns the total number of grid nodes
 **
 ** \return  Total number of nodes
 **
 ** \param[in]   ndim    Space dimension
 ** \param[in]   nxyz    Number of grid nodes along each direction
 **
 *****************************************************************************/
static int st_total_count(int ndim, int *nxyz)
{
  int idim, total;

  total = 1;
  for (idim = 0; idim < ndim; idim++)
    total *= nxyz[idim];
  return (total);
}

/****************************************************************************/
/*!
 **  Calculates the grid extension in a given grid direction
 **
 ** \param[in]   db      Db structure
 ** \param[in]   model   Model structure
 ** \param[in]   simu    ST_FFT structure
 ** \param[in]   percent Percentage of the model variance below which the
 **                      covariance is considered as small enough for dilation
 **
 *****************************************************************************/
static void st_grid_dilate(DbGrid *db, Model *model, ST_FFT *simu, double percent)
{
  double xyz0[3], xyz[3][3];
  int i, j, idx, idy, idz, ndx, ndy, ndz, correct, not_ok, not_ok_dir, ndim,
      indg[3];

  /* Origin of the grid */

  ndim = simu->ndim;
  for (i = 0; i < 3; i++)
    indg[i] = 0;
  grid_to_point(db, indg, nullptr, xyz0);

  /* Location of the elementary end point */

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
      indg[j] = 0;
    indg[i] = 1;
    grid_to_point(db, indg, nullptr, xyz[i]);
  }

  /* Coordinates of the grid vector in the rotated space */

  for (j = 0; j < 3; j++)
    for (i = 0; i < 3; i++)
    {
      xyz[j][i] -= xyz0[i];
      if (j >= ndim) xyz[j][i] = 0.;
    }

  /* Evaluate the count of elementary grid mesh (in each direction) */
  /* for the covariance to become negligeable (<percent)*/

  ndx = ndy = ndz = not_ok = 1;

  while (not_ok)
  {
    not_ok = 0;

    /* Extension along X */

    if (ndim >= 1)
    {
      not_ok_dir = 1;
      while (not_ok_dir)
      {
        correct = 1;
        for (idy = 0; idy < ndy && correct; idy++)
          for (idz = 0; idz < ndz && correct; idz++)
            st_check_correct(model, xyz, ndx, idy, idz, percent, &correct);

        if (correct)
          not_ok_dir = 0;
        else
        {
          ndx++;
          not_ok = 1;
        }
      }
    }

    /* Extension along Y */

    if (ndim >= 2)
    {
      not_ok_dir = 1;
      while (not_ok_dir)
      {
        correct = 1;
        for (idx = 0; idx < ndx && correct; idx++)
          for (idz = 0; idz < ndz && correct; idz++)
            st_check_correct(model, xyz, idx, ndy, idz, percent, &correct);

        if (correct)
          not_ok_dir = 0;
        else
        {
          ndy++;
          not_ok = 1;
        }
      }
    }

    /* Extension along Z */

    if (ndim >= 3)
    {
      not_ok_dir = 1;
      while (not_ok_dir)
      {
        correct = 1;
        for (idx = 0; idx < ndx && correct; idx++)
          for (idy = 0; idy < ndy && correct; idy++)
            st_check_correct(model, xyz, idx, idy, ndz, percent, &correct);

        if (correct)
          not_ok_dir = 0;
        else
        {
          ndz++;
          not_ok = 1;
        }
      }
    }
  }

  /* Ultimate corrections */

  if (ndim < 3) ndz = 0;
  if (ndim < 2) ndy = 0;
  if (ndim < 1) ndx = 0;

  /* Storing arguments */

  simu->shift[0] = ndx;
  simu->shift[1] = ndy;
  simu->shift[2] = ndz;

  /* Optional printout */

  if (OptDbg::query(EDbg::SIMULATE))
  {
    message("Grid Dilation parameters :\n");
    if (ndim >= 1) message("- Number of Nodes along X = %d\n", ndx);
    if (ndim >= 2) message("- Number of Nodes along Y = %d\n", ndy);
    if (ndim >= 3) message("- Number of Nodes along Z = %d\n", ndz);
  }

  return;
}

/****************************************************************************/
/*!
 **  Free core
 **
 ** \param[in]  simu ST_FFT structure
 **
 *****************************************************************************/
static void st_simfft_free(ST_FFT *simu)

{
  simu->cmat = (double*) mem_free((char* ) simu->cmat);
  simu->rnd = (double*) mem_free((char* ) simu->rnd);
  simu->u = (double*) mem_free((char* ) simu->u);
  simu->v = (double*) mem_free((char* ) simu->v);
  simu->sizes_alloc = 0;

  return;
}

/****************************************************************************/
/*!
 **  Dimension the ST_FFT structure
 **
 ** \param[in]  db      Db structure
 ** \param[in]  model   Model structure
 ** \param[in]  percent Percentage of the model variance below which the
 **                     covariance is considered as small enough for dilation
 **
 ** \param[out]  simu   ST_FFT structure
 **
 *****************************************************************************/
static int st_simfft_alloc(DbGrid *db, Model *model, double percent, ST_FFT *simu)
{
  int i, error, ndim, nval;

  /* Initializations */

  error = 1;
  ndim = simu->ndim;

  /* Dilate the grid */

  st_grid_dilate(db, model, simu, percent);

  /* Determine the grid extension */

  for (i = 0; i < 3; i++)
  {
    if (i < ndim)
    {
      nval = st_get_optimal_even_number(simu->shift[i] + db->getNX(i));
      simu->dims[i] = nval;
      simu->dim2[i] = nval / 2;
    }
    else
    {
      simu->dims[i] = 1;
      simu->dim2[i] = 0;
    }
  }
  simu->sizes_alloc = st_total_count(ndim, simu->dims);
  if (OptDbg::query(EDbg::SIMULATE))
  {
    message("Grid parameters after Optimal Dilation :\n");
    if (ndim >= 1) message("- Number of Nodes along X = %d\n", simu->dims[0]);
    if (ndim >= 2) message("- Number of Nodes along Y = %d\n", simu->dims[1]);
    if (ndim >= 3) message("- Number of Nodes along Z = %d\n", simu->dims[2]);
    message("- Total count of grid nodes = %d\n", simu->sizes_alloc);
  }

  /* Core allocation */

  simu->cmat = (double*) mem_alloc(sizeof(double) * simu->sizes_alloc, 0);
  if (simu->cmat == nullptr) goto label_end;
  simu->rnd = (double*) mem_alloc(sizeof(double) * simu->sizes_alloc, 0);
  if (simu->rnd == nullptr) goto label_end;
  simu->u = (double*) mem_alloc(sizeof(double) * simu->sizes_alloc, 0);
  if (simu->u == nullptr) goto label_end;
  simu->v = (double*) mem_alloc(sizeof(double) * simu->sizes_alloc, 0);
  if (simu->v == nullptr) goto label_end;
  for (i = 0; i < simu->sizes_alloc; i++)
  {
    simu->cmat[i] = 0.;
    simu->rnd[i] = 0.;
    simu->u[i] = 0.;
    simu->v[i] = 0.;
  }

  /* Set the return flag */

  error = 0;

  label_end: if (error) st_simfft_free(simu);
  return (error);
}

/****************************************************************************/
/*!
 **  Checks the environment for simulations
 **
 ** \return  Error return code
 **
 ** \param[in]  db     grid Db structure
 ** \param[in]  simu   ST_FFT Structure
 ** \param[in]  model  Model structure
 **
 *****************************************************************************/
static int st_check_simfft_environment(Db *db, ST_FFT *simu, Model *model)
{
  int error = 1;
  int ndim = simu->ndim;

  /**************************************************************/
  /* Check if the Space dimension is compatible with the method */
  /**************************************************************/

  if (ndim < 1 || ndim > 3)
  {
    messerr("The FFT Method is not a relevant simulation model");
    messerr("for this Space Dimension (%d)", ndim);
    return 1;
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
      return 1;
    }
    if (model->getCovaNumber() <= 0)
    {
      messerr("The number of covariance must be positive");
      return 1;
    }
    if (model->getDimensionNumber() <= 0)
    {
      messerr("The Space Dimension must be positive = %d",
              model->getDimensionNumber());
      return 1;
    }
    if (model->getDimensionNumber() != ndim)
    {
      messerr("The Space Dimension of the Db structure (%d)", ndim);
      messerr("Does not correspond to the Space Dimension of the model (%d)",
              model->getDimensionNumber());
      return 1;
    }
  }

  /********************/
  /* Calculate the Db */
  /********************/

  VectorDouble db_mini(ndim);
  VectorDouble db_maxi(ndim);
  db_extension(db, db_mini, db_maxi);

  if (model != nullptr)
    model->setField(ut_vector_extension_diagonal(db_mini, db_maxi));

  /* Set the error return code */

  error = 0;

  label_end:
  return (error);
}

/****************************************************************************/
/*!
 **  Print a plane of values
 **
 ** \param[in]  simu ST_FFT structure
 ** \param[in]  iz   Target plane
 **
 *****************************************************************************/
static void st_print(ST_FFT *simu, int iz)
{
  if (! DEBUG) return;

  print_matrix("U", 0, 1, simu->dims[0], simu->dims[1], NULL, &U(0, 0, iz));
  print_matrix("V", 0, 1, simu->dims[0], simu->dims[1], NULL, &V(0, 0, iz));
}

/****************************************************************************/
/*!
 **  Operate the symmetry for a 1-D space
 **
 ** \param[in]  simu ST_FFT structure
 **
 *****************************************************************************/
static void st_simfft_sym_1(ST_FFT *simu)

{
  int ix, jx;

  // A(1) and A(N1/2+1) are real
  for (ix = 0; ix < simu->dims[0]; ix += simu->dim2[0])
    st_set_zero(simu, ix, 0, 0);

  // A(j) = A*(N1-j+2) for j in [2,N1/2]
  for (ix = 1; ix < simu->dim2[0]; ix++)
  {
    jx = simu->dims[0] - ix;
    st_set_conj(simu, ix, 0, 0, jx, 0, 0);
  }

  st_print(simu, 0);
  return;
}

/****************************************************************************/
/*!
 **  Operate the symmetry for a 2-D space
 **
 ** \param[in]  simu  ST_FFT structure
 ** \param[in]  iz0   fixed third index
 **
 *****************************************************************************/
static void st_simfft_sym_2(ST_FFT *simu, int iz0)
{
  int ix, iy, jx, jy;

  // A(1,1), A(N1/2,1), A(1,N2/2) and A(N1/2,N2/2) real
  for (iy = 0; iy < simu->dims[1]; iy += simu->dim2[1])
    for (ix = 0; ix < simu->dims[0]; ix += simu->dim2[0])
      st_set_zero(simu, ix, iy, iz0);

  // A(1,k)      = A*(1,N2-k+2)      for k in [2,N2/2]
  // A(N2/2+1,k) = A*(N2/2+1,N2-k+2) for k in [2,N2/2]
  for (iy = 1; iy < simu->dim2[1]; iy++)
    for (ix = 0; ix < simu->dims[0]; ix += simu->dim2[0])
    {
      jy = simu->dims[1] - iy;
      st_set_conj(simu, ix, iy, iz0, ix, jy, iz0);
    }

  // A(j,1)      = A*(N1-j+2,1)      for j in [2,N1/2]
  // A(j,N2/2+1) = A*(N1-j+2,N2/2+1) for j in [2,N1/2]
  for (iy = 0; iy < simu->dims[1]; iy += simu->dim2[1])
    for (ix = 1; ix < simu->dim2[0]; ix++)
    {
      jx = simu->dims[0] - ix;
      st_set_conj(simu, ix, iy, iz0, jx, iy, iz0);
    }

  // A(j,k) = A*(N1-j+2,N2-k+2) for j in [2,N1/2] and k in [2,N2/2]
  for (iy = 1; iy < simu->dim2[1]; iy++)
    for (ix = 1; ix < simu->dim2[0]; ix++)
    {
      jx = simu->dims[0] - ix;
      jy = simu->dims[1] - iy;
      st_set_conj(simu, ix, iy, iz0, jx, jy, iz0);
    }

  // A(j,N2-k+2) = A*(N1-j+2,k) for j in [2,N1/2] and k in [2,N2/2]
  for (iy = 1; iy < simu->dim2[1]; iy++)
    for (ix = 1; ix < simu->dim2[0]; ix++)
    {
      jx = simu->dims[0] - ix;
      jy = simu->dims[1] - iy;
      st_set_conj(simu, ix, jy, iz0, jx, iy, iz0);
    }

  st_print(simu, 0);
  return;
}

/****************************************************************************/
/*!
 **  Operate the symmetry for a 3-D space
 **
 ** \param[in]  simu  ST_FFT structure
 **
 *****************************************************************************/
static void st_simfft_sym_3(ST_FFT *simu)

{
  int ix, iy, iz, jx, jy, jz;

  // For l=1 or N3/2+1, use the 2-D symmetry
  for (iz = 0; iz < simu->dims[2]; iz += simu->dim2[2])
    st_simfft_sym_2(simu, iz);

  // For the other planes:

  // A(1,1,l)           = A*(1,1,N3-l+2)           for l in [2,N3/2]
  // A(N1/2+1,1,l)      = A*(N1/2+1,1,N3-l+2)      for l in [2,N3/2]
  // A(1,N2/2+1,l)      = A*(1,N2/2+1,N3-l+2)      for l in [2,N3/2]
  // A(N1/2+1,N2/2+1,l) = A*(N1/2+1,N2/2+1,N3-l+2) for l in [2,N3/2]
  for (iz = 1; iz < simu->dim2[2]; iz++)
    for (iy = 0; iy < simu->dims[1]; iy += simu->dim2[1])
      for (ix = 0; ix < simu->dims[0]; ix += simu->dim2[0])
      {
        jz = simu->dims[2] - iz;
        st_set_conj(simu, ix, iy, iz, ix, iy, jz);
      }

  // A(1,k,l)      = A*(1,N2-k+2,N3-l+2)      for k in [2,N2/2] and l in [2,N3/2]
  // A(N1/2+1,k,l) = A*(N1/2+1,N2-k+2,N3-l+2) for k in [2,N2/2] and l in [2,N3/2]
  for (iz = 1; iz < simu->dim2[2]; iz++)
    for (iy = 1; iy < simu->dim2[1]; iy++)
      for (ix = 0; ix < simu->dims[0]; ix += simu->dim2[0])
      {
        jy = simu->dims[1] - iy;
        jz = simu->dims[2] - iz;
        st_set_conj(simu, ix, iy, iz, ix, jy, jz);
      }

  // A(1,N2-k+2,l)      = A*(1,k,N3/-l+2)      for k in [2,N2/2] and l in [2,N3/2]
  // A(N1/2+1,N2-k+2,l) = A*(N1/2+1,k,N3/-l+2) for k in [2,N2/2] and l in [2,N3/2]
  for (iz = 1; iz < simu->dim2[2]; iz++)
    for (iy = 1; iy < simu->dim2[1]; iy++)
      for (ix = 0; ix < simu->dims[0]; ix += simu->dim2[0])
      {
        jy = simu->dims[1] - iy;
        jz = simu->dims[2] - iz;
        st_set_conj(simu, ix, jy, iz, ix, iy, jz);
      }

  // A(j,1,l)      = A*(N1-j+2,1,N3-l+2)      for j in [2,N1/2] and l in [2,N3/2]
  // A(j,N2/2+1,l) = A*(N1-j+2,N2/2+1,N3-l+2) for j in [2,N1/2] and l in [2,N3/2]
  for (iz = 1; iz < simu->dim2[2]; iz++)
    for (iy = 0; iy < simu->dims[1]; iy += simu->dim2[1])
      for (ix = 1; ix < simu->dim2[0]; ix++)
      {
        jx = simu->dims[0] - ix;
        jz = simu->dims[2] - iz;
        st_set_conj(simu, ix, iy, iz, jx, iy, jz);
      }

  // A(N1-j+2,1,l)      = A*(j,1,N3-l+2)      for j in [2,N1/2] and l in [2,N3/2]
  // A(N1-j+2,N2/2+1,l) = A*(j,N2/2+1,N3-l+2) for j in [2,N1/2] and l in [2,N3/2]
  for (iz = 1; iz < simu->dim2[2]; iz++)
    for (iy = 0; iy < simu->dims[1]; iy += simu->dim2[1])
      for (ix = 1; ix < simu->dim2[0]; ix++)
      {
        jx = simu->dims[0] - ix;
        jz = simu->dims[2] - iz;
        st_set_conj(simu, jx, iy, iz, ix, iy, jz);
      }

  // A(j,k,l) = A*(N1-j+2,N2-k+2,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
  for (iz = 1; iz < simu->dim2[2]; iz++)
    for (iy = 1; iy < simu->dim2[1]; iy++)
      for (ix = 1; ix < simu->dim2[0]; ix++)
      {
        jx = simu->dims[0] - ix;
        jy = simu->dims[1] - iy;
        jz = simu->dims[2] - iz;
        st_set_conj(simu, ix, iy, iz, jx, jy, jz);
      }

  // A(N1-j+2,N2-k+2,l) = A*(j,k,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
  for (iz = 1; iz < simu->dim2[2]; iz++)
    for (iy = 1; iy < simu->dim2[1]; iy++)
      for (ix = 1; ix < simu->dim2[0]; ix++)
      {
        jx = simu->dims[0] - ix;
        jy = simu->dims[1] - iy;
        jz = simu->dims[2] - iz;
        st_set_conj(simu, jx, jy, iz, ix, iy, jz);
      }

  // A(j,N2-k+2,l) = A*(N2-j+2,k,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
  for (iz = 1; iz < simu->dim2[2]; iz++)
    for (iy = 1; iy < simu->dim2[1]; iy++)
      for (ix = 1; ix < simu->dim2[0]; ix++)
      {
        jx = simu->dims[0] - ix;
        jy = simu->dims[1] - iy;
        jz = simu->dims[2] - iz;
        st_set_conj(simu, ix, jy, iz, jx, iy, jz);
      }

  // A(N1-j+2,k,l) = A*(j,k,N3-l+2) for j in [2,N1/2], k in [2,N2/2] and l in [2,N3/2]
  for (iz = 1; iz < simu->dim2[2]; iz++)
    for (iy = 1; iy < simu->dim2[1]; iy++)
      for (ix = 1; ix < simu->dim2[0]; ix++)
      {
        jx = simu->dims[0] - ix;
        jy = simu->dims[1] - iy;
        jz = simu->dims[2] - iz;
        st_set_conj(simu, jx, iy, iz, ix, jy, jz);
      }

  for (iz = 0; iz < simu->dims[2]; iz++)
    st_print(simu, iz);
  return;
}

/****************************************************************************/
/*!
 **  Operate the symmetry
 **
 ** \param[in]  simu  ST_FFT structure
 **
 *****************************************************************************/
static void st_simfft_symmetry(ST_FFT *simu)
{

  /* Dispatch according to the space dimension */

  switch (simu->ndim)
  {
    case 1:
      st_simfft_sym_1(simu);
      break;

    case 2:
      st_simfft_sym_2(simu, 0);
      break;

    case 3:
      st_simfft_sym_3(simu);
      break;
  }

  return;
}

/****************************************************************************/
/*!
 **  Perform a non-conditional simulation on the grid
 **
 ** \param[in]  db    Db structure
 ** \param[in]  simu  ST_FFT structure
 ** \param[in]  iad   address for writing the simulation
 **
 *****************************************************************************/
static void st_simfft_final(DbGrid *db, ST_FFT *simu, int iad)
{
  int ix, iy, iz, jx, jy, jz, ecr;

  /* Perform the Inverse Fast Fourier Transform */

  (void) fftn(simu->ndim, simu->dims, simu->u, simu->v, 1, 1.);

  /* Retrieving the simulation */

  for (iz = ecr = 0; iz < db->getNX(2); iz++)
    for (iy = 0; iy < db->getNX(1); iy++)
      for (ix = 0; ix < db->getNX(0); ix++, ecr++)
      {
        jx = ix + simu->shift[0];
        jy = iy + simu->shift[1];
        jz = iz + simu->shift[2];
        db->updArray(ecr, iad, 4, U(jx, jy, jz));
      }

  return;
}

/****************************************************************************/
/*!
 **  Perform the non-conditional simulation by FFT method on a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  model   Model structure
 ** \param[in]  seed    Value of the seed
 ** \param[in]  nbsimu  Number of simulations
 ** \param[in]  percent Percentage of the model variance below which the
 **                      covariance is considered as small enough for dilation
 ** \param[in]  flag_aliasing  1 for anti-aliasing procedure; 0 otherwise
 **
 *****************************************************************************/
int simfft_f(DbGrid *db,
             Model *model,
             int seed,
             int nbsimu,
             double percent,
             int flag_aliasing)
{
  ST_FFT simu;
  int isimu, iptr, error;

  /* Initializations */

  error = 1;
  FLAG_ALIASING = flag_aliasing;
  if (seed != 0) law_set_random_seed(seed);
  st_simfft_init(db, &simu);
  if (st_check_simfft_environment(db, &simu, model)) goto label_end;

  /* Add the attributes for storing the results in the data base */

  iptr = db->addColumnsByConstant(nbsimu, 0.);

  /* Construction of the Simu_FFT structure and core allocation */

  if (st_simfft_alloc(db, model, percent, &simu)) goto label_end;

  /* Preparation of the FFT environment */

  (void) st_simfft_prepar(db, model, &simu, 1);

  /* Processing */

  for (isimu = 0; isimu < nbsimu; isimu++)
  {

    /* Initiate the random normal values */

    st_simfft_random(&simu);

    /* Apply the symmetry */

    st_simfft_symmetry(&simu);

    /* Perform the simulation */

    st_simfft_final(db, &simu, iptr + isimu);
  }

  /* Set the error code */

  error = 0;

  label_end: st_simfft_free(&simu);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the exponential of the scaled correlation
 **
 ** \return  Value of the transformed correlation
 **
 ** \param[in]  simu   ST_FFT structure
 ** \param[in]  sigma  Logarithmic variance value
 ** \param[in]  ix     Index for the discretized covariance along X
 ** \param[in]  iy     Index for the discretized covariance along Y
 ** \param[in]  iz     Index for the discretized covariance along Z
 **
 *****************************************************************************/
static double st_rho_sigma(ST_FFT *simu, double sigma, int ix, int iy, int iz)
{
  double rho;

  rho = simu->cmat[IND(ix, iy, iz)];
  if (!FFFF(sigma)) rho = exp(sigma * sigma * rho);
  return (rho);
}

/****************************************************************************/
/*!
 **  Calculate the mean lognormal covariance over the block (1-D)
 **
 ** \return  Mean lognormal covariance
 **
 ** \param[in]  simu   ST_FFT structure
 ** \param[in]  sigma  Logarithmic variance value
 **
 *****************************************************************************/
static double st_support_1(ST_FFT *simu, double sigma)
{
  int ix, iix;
  double value, rho;

  value = 0.;
  for (ix = -simu->nx[0]; ix <= simu->nx[0]; ix++)
  {
    iix = (ix < 0) ? simu->dims[0] + ix :
                     ix;
    rho = st_rho_sigma(simu, sigma, iix, 0, 0);
    value += (simu->nx[0] - ABS(ix)) * rho;
  }
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the mean lognormal covariance over the block (2-D)
 **
 ** \return  Mean lognormal covariance
 **
 ** \param[in]  simu   ST_FFT structure
 ** \param[in]  sigma  Logarithmic variance value
 **
 *****************************************************************************/
static double st_support_2(ST_FFT *simu, double sigma)
{
  int ix, iy, iix, iiy;
  double value, rho;

  value = 0.;
  for (ix = -simu->nx[0]; ix <= simu->nx[0]; ix++)
    for (iy = -simu->nx[1]; iy <= simu->nx[1]; iy++)
    {
      iix = (ix < 0) ? simu->dims[0] + ix :
                       ix;
      iiy = (iy < 0) ? simu->dims[1] + iy :
                       iy;
      rho = st_rho_sigma(simu, sigma, iix, iiy, 0);
      value += ((simu->nx[0] - ABS(ix)) * (simu->nx[1] - ABS(iy)) * rho);
    }
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the mean lognormal covariance over the block (3-D)
 **
 ** \return  Mean lognormal covariance
 **
 ** \param[in]  simu   ST_FFT structure
 ** \param[in]  sigma  Logarithmic variance value
 **
 *****************************************************************************/
static double st_support_3(ST_FFT *simu, double sigma)
{
  int ix, iy, iz, iix, iiy, iiz;
  double value, rho;

  value = 0.;
  for (ix = -simu->nx[0]; ix <= simu->nx[0]; ix++)
    for (iy = -simu->nx[1]; iy <= simu->nx[1]; iy++)
      for (iz = -simu->nx[2]; iz <= simu->nx[2]; iz++)
      {
        iix = (ix < 0) ? simu->dims[0] + ix :
                         ix;
        iiy = (iy < 0) ? simu->dims[1] + iy :
                         iy;
        iiz = (iz < 0) ? simu->dims[2] + iz :
                         iz;
        rho = st_rho_sigma(simu, sigma, iix, iiy, iiz);
        value += ((simu->nx[0] - ABS(ix)) * (simu->nx[1] - ABS(iy))
                  * (simu->nx[2] - ABS(iz)) * rho);
      }
  return (value);
}

/****************************************************************************/
/*!
 **  Calculate the mean lognormal covariance over the block
 **
 ** \return  Mean lognormal covariance
 **
 ** \param[in]  simu   ST_FFT structure
 ** \param[in]  sigma  Logarithmic variance value
 **
 *****************************************************************************/
static double st_support(ST_FFT *simu, double sigma)
{
  double value, scale;
  int idim;

  value = 0.;
  if (sigma == 0.) return (TEST);

  switch (simu->ndim)
  {
    case 1:
      value = st_support_1(simu, sigma);
      break;

    case 2:
      value = st_support_2(simu, sigma);
      break;

    case 3:
      value = st_support_3(simu, sigma);
      break;
  }

  /* Calculate the scale */

  scale = 1.;
  for (idim = 0; idim < simu->ndim; idim++)
    scale *= (simu->nx[idim] * simu->nx[idim]);
  value /= scale;

  /* Back-transform into change of support coefficient */

  if (!FFFF(sigma)) value = log(value) / (sigma * sigma);

  return (sqrt(value));
}

/****************************************************************************/
/*!
 **  Calculate the change of support coefficients by FFT method
 **  in the lognormal case on a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  model   Model structure
 ** \param[in]  percent Percentage of the model variance below which the
 **                      covariance is considered as small enough for dilation
 ** \param[in]  flag_aliasing  1 for anti-aliasing procedure; 0 otherwise
 ** \param[in]  nval    Number of logarithmic variances
 ** \param[in]  sigma   Array of logarithmic variances
 **
 ** \param[out] r2val   r^2 coefficients
 ** \param[out] coeffs  r^2 coefficients for given logarithmic variances
 **
 *****************************************************************************/
int simfft_support(DbGrid *db,
                   Model *model,
                   double percent,
                   int flag_aliasing,
                   int nval,
                   double *sigma,
                   double *r2val,
                   double *coeffs)
{
  ST_FFT simu;
  int error, ival;

  /* Initializations */

  error = 1;
  FLAG_ALIASING = flag_aliasing;
  st_simfft_init(db, &simu);
  if (st_check_simfft_environment(db, &simu, model)) goto label_end;

  /* Construction of the Simu_FFT structure and core allocation */

  if (st_simfft_alloc(db, model, percent, &simu)) goto label_end;

  /* Preparation of the FFT environment */

  (void) st_simfft_prepar(db, model, &simu, 0);

  /* Calculate the correlation matrix (possibly rescaled) */

  (void) fftn(simu.ndim, simu.dims, simu.cmat, simu.rnd, 1, 1.);
  if (DEBUG)
    print_matrix("Discretized Covariance (after scaling)", 0, 1, simu.dims[0],
                 simu.dims[1], NULL, simu.cmat);

  /* Loop on the different lognormal variances */

  *r2val = st_support(&simu, TEST);
  for (ival = 0; ival < nval; ival++)
    coeffs[ival] = st_support(&simu, sigma[ival]);

  /* Set the error code */

  error = 0;

  label_end: st_simfft_free(&simu);
  return (error);
}
