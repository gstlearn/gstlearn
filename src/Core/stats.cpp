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
#include "geoslib_old_f.h"

#include "Morpho/Morpho.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Stats/Classical.hpp"

#include <math.h>
#include <string.h>

/*! \cond */
#define G_ADDRESS(ix,iy,iz,nxyz)    ((ix) + nxyz[0] * ((iy) + nxyz[1] * (iz)))
#define N1_TAB(ix,iy,iz) (numtab1[G_ADDRESS(ix,iy,iz,nxyz1)])
#define N2_TAB(ix,iy,iz) (numtab2[G_ADDRESS(ix,iy,iz,nxyz2)])
#define V1_TAB(ix,iy,iz) (valtab1[G_ADDRESS(ix,iy,iz,nxyz1)])
#define V2_TAB(ix,iy,iz) (valtab2[G_ADDRESS(ix,iy,iz,nxyz2)])
#define D1_TAB(ix,iy,iz) (N1_TAB(ix,iy,iz) > 0 &&     \
                          ! FFFF(V1_TAB(ix,iy,iz)) && \
                          V1_TAB(ix,iy,iz) > 0)
#define RESIDUALS(icut,iech) (residuals[(icut) * nech + (iech)])
#define NBGH(ivois,idim)     (nbgh[ndim * (ivois) + (idim)])
#define TABINI(iseed,idim)   (tabini[ndim * (iseed) + (idim)])
#define TABCUR(iseed,idim)   (tabcur[ndim * (iseed) + (idim)])
#define TRAJEC(iseed,iter,idim) (trsave[(niter * (iseed) + (iter)) * ndim + (idim)])
/*! \endcond */

static int DEBUG = 0;

/****************************************************************************/
/*!
 **  Load the subgrid from the Input Db
 **
 ** \return Return the total probability of finding joins
 **
 ** \param[in]  verbose   Verbose flag
 ** \param[in]  flag_ffff 1 replace masked values by 0
 **                       0 replace masked values by FFFF
 ** \param[in]  iech0     Rank of the output grid cell
 ** \param[in]  nech0     Number of cells of the output grid
 ** \param[in]  ntot      Dimension of arrays 'numtab1' and 'valtab1'
 ** \param[in]  dbgrid    Db for the input grid
 ** \param[in]  ind0      Origin of the Output Db within the Input Db
 ** \param[in]  nxyz      Mesh of the Output Db (expressed in Input Db)
 ** \param[in]  ixyz      Indices of the starting node of the Output Db
 **
 ** \param[out] numtab1   Array containing the sample count
 ** \param[out] valtab1   Array containing the sample value
 **
 ** \remarks  The array ind0, ixyz and nxyz are dimensioned to ndim
 **
 *****************************************************************************/
static double st_extract_subgrid(int verbose,
                                 int flag_ffff,
                                 int iech0,
                                 int nech0,
                                 int ntot,
                                 DbGrid *dbgrid,
                                 int *ind0,
                                 int *ixyz,
                                 int *nxyz,
                                 double *numtab1,
                                 double *valtab1)
{
  int ix, iy, iz, jx, jy, jz, ind, ecr, ndim;
  double proba, value;

  /* Initializations */

  ndim = dbgrid->getNDim();
  VectorInt iwork2(ndim);
  for (int i = 0; i < ntot; i++)
  {
    numtab1[i] = 0;
    valtab1[i] = 0.;
  }

  for (int idim = 0; idim < 3; idim++)
  {
    if (idim < ndim) continue;
    ixyz[idim] = 0;
    nxyz[idim] = 1;
    ind0[idim] = 0;
  }

  ecr = 0;
  proba = 0.;
  for (iz = 0; iz < nxyz[2]; iz++)
    for (iy = 0; iy < nxyz[1]; iy++)
      for (ix = 0; ix < nxyz[0]; ix++)
      {

        /* Get the address of a sample of the subgrid */

        jx = ind0[0] + ixyz[0] * nxyz[0] + ix;
        if (jx < 0 || jx > dbgrid->getNX(0)) continue;
        jy = ind0[1] + ixyz[1] * nxyz[1] + iy;
        if (jy < 0 || jy > dbgrid->getNX(1)) continue;
        jz = ind0[2] + ixyz[2] * nxyz[2] + iz;
        if (jz < 0 || jz > dbgrid->getNX(2)) continue;

        /* Get the node index within the Input Db */

        if (ndim >= 1) iwork2[0] = jx;
        if (ndim >= 2) iwork2[1] = jy;
        if (ndim >= 3) iwork2[2] = jz;
        ind = db_index_grid_to_sample(dbgrid, iwork2.data());
        numtab1[ecr] = 1.;
        value = dbgrid->isActive(ind) ? dbgrid->getLocVariable(ELoc::Z,ind, 0) :
                                        TEST;
        if (FFFF(value))
          valtab1[ecr] = (flag_ffff) ? 0 :
                                       TEST;
        else
        {
          valtab1[ecr] = value;
          proba += value;
        }
        ecr++;
      }

  /* Optional verbose option */

  if (verbose)
  {
    message("Output cell %3d/%3d = %d", iech0 + 1, nech0, nxyz[0]);
    for (int idim = 1; idim < ndim; idim++)
      message("x%d", nxyz[idim]);
    message(" cells of Input Grid (Proba=%lf)\n", proba);
  }
  return (proba);
}

/****************************************************************************/
/*!
 **  Divide by 2 in integer (upper rounded value)
 **
 ** \return  1 if the input value is larger than 2; 0 otherwise
 **
 ** \param[in]  nxyz      Dimensions of the subgrid
 ** \param[in]  orient    Rank of the target direction
 **
 *****************************************************************************/
static int st_divide_by_2(int *nxyz, int orient)
{
  int ival;

  ival = nxyz[orient];
  if (ival <= 1) return (0);

  ival = (int) floor((double) (ival + 1.) / 2);
  nxyz[orient] = ival;
  return (1);
}

/****************************************************************************/
/*!
 **  Perform the arithmetic mean
 **
 ** \param[in]  idim      Direction of calculation
 ** \param[in]  nxyz1     Dimensions of the initial subgrid
 ** \param[in]  nxyz2     Dimensions of the final subgrid
 **
 ** \param[out] numtab1   Array containing the sample count
 ** \param[out] numtab2   Array containing the sample count
 ** \param[out] valtab1   Array containing the sample value
 ** \param[out] valtab2   Array containing the sample value
 **
 *****************************************************************************/
static void st_mean_arith(int idim,
                          const int *nxyz1,
                          const int *nxyz2,
                          const double *numtab1,
                          double *numtab2,
                          double *valtab1,
                          double *valtab2)
{
  int ix, iy, iz, ix1, ix2, iy1, iy2, iz1, iz2;

  for (iz = 0; iz < nxyz2[2]; iz++)
    for (iy = 0; iy < nxyz2[1]; iy++)
      for (ix = 0; ix < nxyz2[0]; ix++)
      {
        N2_TAB(ix,iy,iz) = V2_TAB(ix,iy,iz) = 0.;
        switch (idim)
        {
          case 0:
            ix1 = 2 * ix;
            ix2 = 2 * ix + 1;
            if (D1_TAB(ix1, iy, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix1, iy, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix1,iy,iz) * V1_TAB(ix1, iy, iz);
            }
            if (ix2 < nxyz1[0] && D1_TAB(ix2, iy, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix2, iy, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix2,iy,iz) * V1_TAB(ix2, iy, iz);
            }
            break;

          case 1:
            iy1 = 2 * iy;
            iy2 = 2 * iy + 1;
            if (D1_TAB(ix, iy1, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy1, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy1,iz) * V1_TAB(ix, iy1, iz);
            }
            if (iy2 < nxyz1[1] && D1_TAB(ix, iy2, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy2, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy2,iz) * V1_TAB(ix, iy2, iz);
            }
            break;

          case 2:
            iz1 = 2 * iz;
            iz2 = 2 * iz + 1;
            if (D1_TAB(ix, iy, iz1))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy, iz1);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy,iz1) * V1_TAB(ix, iy, iz1);
            }
            if (iz2 < nxyz1[2] && D1_TAB(ix, iy, iz2))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy, iz2);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy,iz2) * V1_TAB(ix, iy, iz2);
            }
            break;
        }
        V2_TAB(ix,iy,iz) =
            (N2_TAB(ix,iy,iz) > 0) ?
                                     V2_TAB(ix,iy,iz) / N2_TAB(ix, iy, iz) :
                                     TEST;
      }
}

/****************************************************************************/
/*!
 **  Perform the harmonic mean
 **
 ** \param[in]  idim      Direction of calculation
 ** \param[in]  nxyz1     Dimensions of the initial subgrid
 ** \param[in]  nxyz2     Dimensions of the final subgrid
 **
 ** \param[out] numtab1   Array containing the input  sample count
 ** \param[out] numtab2   Array containing the output sample count
 ** \param[out] valtab1   Array containing the input  sample value
 ** \param[out] valtab2   Array containing the output sample value
 **
 *****************************************************************************/
static void st_mean_harmo(int idim,
                          const int *nxyz1,
                          const int *nxyz2,
                          const double *numtab1,
                          double *numtab2,
                          double *valtab1,
                          double *valtab2)
{
  int ix, iy, iz, ix1, ix2, iy1, iy2, iz1, iz2;

  for (iz = 0; iz < nxyz2[2]; iz++)
    for (iy = 0; iy < nxyz2[1]; iy++)
      for (ix = 0; ix < nxyz2[0]; ix++)
      {
        N2_TAB(ix,iy,iz) = V2_TAB(ix,iy,iz) = 0.;
        switch (idim)
        {
          case 0:
            ix1 = 2 * ix;
            ix2 = 2 * ix + 1;
            if (D1_TAB(ix1, iy, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix1, iy, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix1,iy,iz) / V1_TAB(ix1, iy, iz);
            }
            if (ix2 < nxyz1[0] && D1_TAB(ix2, iy, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix2, iy, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix2,iy,iz) / V1_TAB(ix2, iy, iz);
            }
            break;

          case 1:
            iy1 = 2 * iy;
            iy2 = 2 * iy + 1;
            if (D1_TAB(ix, iy1, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy1, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy1,iz) / V1_TAB(ix, iy1, iz);
            }
            if (iy2 < nxyz1[1] && D1_TAB(ix, iy2, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy2, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy2,iz) / V1_TAB(ix, iy2, iz);
            }
            break;

          case 2:
            iz1 = 2 * iz;
            iz2 = 2 * iz + 1;
            if (D1_TAB(ix, iy, iz1))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy, iz1);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy,iz1) / V1_TAB(ix, iy, iz1);
            }
            if (iz2 < nxyz1[2] && D1_TAB(ix, iy, iz2))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy, iz2);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy,iz2) / V1_TAB(ix, iy, iz2);
            }
            break;
        }
        V2_TAB(ix, iy, iz) = (ABS(V2_TAB(ix, iy, iz)) > 1.e-10)
                               ? N2_TAB(ix, iy, iz) / V2_TAB(ix, iy, iz)
                               : TEST;
      }
}

/****************************************************************************/
/*!
 **  Update the dimensions and calculate the current number of cells
 **
 ** \return Number of cells
 **
 ** \param[in]  nxyz1     Dimensions of the initial subgrid
 ** \param[in]  numtab1   Array containing the sample count
 ** \param[in]  valtab1   Array containing the sample value
 **
 ** \param[out] nxyz2     Dimensions of the final subgrid
 ** \param[out] numtab2   Array containing the sample count
 ** \param[out] valtab2   Array containing the sample value
 **
 *****************************************************************************/
static int st_recopy(const int *nxyz1,
                     const double *numtab1,
                     const double *valtab1,
                     int *nxyz2,
                     double *numtab2,
                     double *valtab2)
{
  int i, ncell;

  /* Update the dimension */

  ncell = 1;
  for (i = 0; i < 3; i++)
  {
    nxyz2[i] = nxyz1[i];
    ncell *= nxyz1[i];
  }

  /* Update the contents of the arrays */

  for (i = 0; i < ncell; i++)
  {
    numtab2[i] = numtab1[i];
    valtab2[i] = valtab1[i];
  }

  return (ncell);
}

/****************************************************************************/
/*!
 **  Print the generated grids (for counts and values)
 **
 ** \param[in]  subtitle  Subtitle
 ** \param[in]  nxyz      Dimensions of the grid
 **
 ** \param[out] numtab    Array containing the sample count
 ** \param[out] valtab    Array containing the sample value
 **
 ****************************************************************************/
static void st_print_grid(const char *subtitle,
                          int nxyz[3],
                          double *numtab,
                          double *valtab)
{
  char string[100];
  int iz, shift;

  /* Initializations */

  shift = nxyz[0] * nxyz[1];

  /* Loop on the third dimension */

  for (iz = 0; iz < nxyz[2]; iz++)
  {
    (void) gslSPrintf(string, "%s Values (iz=%d)\n", subtitle, iz + 1);
    message(string);
    print_matrix(NULL, 0, 0, nxyz[0], nxyz[1], NULL, &valtab[iz * shift]);
    (void) gslSPrintf(string, "%s Counts (iz=%d)\n", subtitle, iz + 1);
    message(string);
    print_matrix(NULL, 0, 0, nxyz[0], nxyz[1], NULL, &numtab[iz * shift]);
  }
  message("\n");
}

/****************************************************************************/
/*!
 **  Print statistics after basic upscaling operation
 **
 ** \param[in]  title     Title for the printout
 ** \param[in]  nxyz      Grid dimension
 ** \param[out] valtab    Array containing the sample value
 **
 *****************************************************************************/
static void st_print_upscale(const char *title, int *nxyz, const double *valtab)
{
  double mini, maxi, value;
  int lec, ndef;

  mini = 1.e30;
  maxi = -1.e30;
  lec = ndef = 0;
  for (int iz = 0; iz < nxyz[2]; iz++)
    for (int iy = 0; iy < nxyz[1]; iy++)
      for (int ix = 0; ix < nxyz[0]; ix++)
      {
        value = valtab[lec++];
        if (FFFF(value)) continue;
        if (value < mini) mini = value;
        if (value > maxi) maxi = value;
        ndef++;
      }

  message("%11s ", title);
  message("(%3d %3d %3d) : ", nxyz[0], nxyz[1], nxyz[2]);
  message("%lf %lf (Def=%d)\n", mini, maxi, ndef);
}

/****************************************************************************/
/*!
 **  Upscale the variable defined on a subgrid into one value
 **
 ** \param[in]  orient    Upscaling orientation (0 to 2)
 ** \param[in]  nxyz      Dimensions of the subgrid
 ** \param[in]  flag_save If the array must be saved for keypair
 **
 ** \param[out] numtab0   Array containing the sample count
 ** \param[out] numtab1   Array containing the sample count
 ** \param[out] numtab2   Array containing the sample count
 ** \param[out] valtab0   Array containing the sample value
 ** \param[out] valtab1   Array containing the sample value
 ** \param[out] valtab2   Array containing the sample value
 ** \param[out] res1      First result (Harmonic first)
 ** \param[out] res2      First result (Arithmetic first)
 **
 *****************************************************************************/
static void st_upscale(int orient,
                       int *nxyz,
                       int flag_save,
                       double *numtab0,
                       double *numtab1,
                       double *numtab2,
                       double *valtab0,
                       double *valtab1,
                       double *valtab2,
                       double *res1,
                       double *res2)
{
  int idim, nxyz1[3], nxyz2[3], ncell, flag_debug;

  /* Initializations */

  flag_debug = OptDbg::query(EDbg::UPSCALE);

  /**************************************/
  /* Getting the minimum upscaled value */
  /**************************************/

  if (flag_debug)
  {
    mestitle(1, "Looking for the Minimum Upscaled Value");
    st_print_grid("Initial", nxyz, numtab0, valtab0);
  }
  ncell = st_recopy(nxyz, numtab0, valtab0, nxyz1, numtab1, valtab1);
  ncell = st_recopy(nxyz, numtab0, valtab0, nxyz2, numtab2, valtab2);
  if (flag_save) st_print_upscale("Initial", nxyz1, valtab1);

  while (ncell > 1)
  {

    /* Harmonic mean in the flow direction */

    if (st_divide_by_2(nxyz2, orient))
    {
      st_mean_harmo(orient, nxyz1, nxyz2, numtab1, numtab2, valtab1, valtab2);
      ncell = st_recopy(nxyz2, numtab2, valtab2, nxyz1, numtab1, valtab1);
      if (flag_save) st_print_upscale("Harmonic", nxyz1, valtab1);
      if (flag_debug) st_print_grid("Harmonic", nxyz1, numtab1, valtab1);
    }

    /* Arithmetic mean orthogonal to the flow direction */

    for (idim = 0; idim < 3; idim++)
    {
      if (idim == orient) continue;
      if (st_divide_by_2(nxyz2, idim))
      {
        st_mean_arith(idim, nxyz1, nxyz2, numtab1, numtab2, valtab1, valtab2);
        ncell = st_recopy(nxyz2, numtab2, valtab2, nxyz1, numtab1, valtab1);
        if (flag_save) st_print_upscale("Arithmetic", nxyz1, valtab1);
        if (flag_debug) st_print_grid("Arithmetic", nxyz1, numtab1, valtab1);
      }
    }
  }
  *res1 = valtab1[0];

  /**************************************/
  /* Getting the maximum upscaled value */
  /**************************************/

  if (flag_debug)
  {
    mestitle(1, "Looking for the Maximum Upscaled Value\n");
    st_print_grid("Initial", nxyz, numtab0, valtab0);
  }
  ncell = st_recopy(nxyz, numtab0, valtab0, nxyz1, numtab1, valtab1);
  ncell = st_recopy(nxyz, numtab0, valtab0, nxyz2, numtab2, valtab2);
  if (flag_save) st_print_upscale("Initial", nxyz1, valtab1);

  while (ncell > 1)
  {

    /* Arithmetic mean orthogonal to the flow direction */

    for (idim = 0; idim < 3; idim++)
    {
      if (idim == orient) continue;
      if (st_divide_by_2(nxyz2, idim))
      {
        st_mean_arith(idim, nxyz1, nxyz2, numtab1, numtab2, valtab1, valtab2);
        ncell = st_recopy(nxyz2, numtab2, valtab2, nxyz1, numtab1, valtab1);
        if (flag_save) st_print_upscale("Arithmetic", nxyz1, valtab1);
        if (flag_debug) st_print_grid("Arithmetic", nxyz1, numtab1, valtab1);
      }
    }

    /* Harmonic mean in the flow direction */

    if (st_divide_by_2(nxyz2, orient))
    {
      st_mean_harmo(orient, nxyz1, nxyz2, numtab1, numtab2, valtab1, valtab2);
      ncell = st_recopy(nxyz2, numtab2, valtab2, nxyz1, numtab1, valtab1);
      if (flag_save) st_print_upscale("Harmonic", nxyz1, valtab1);
      if (flag_debug) st_print_grid("Harmonic", nxyz1, numtab1, valtab1);
    }
  }
  *res2 = valtab1[0];
}

/****************************************************************************/
/*!
 **  Check if the first grid is a subgrid of the second one
 **
 ** \return  1 if grids are multiple
 **
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  title      Title of the work (used only for verbose case)
 ** \param[in]  dbgrid1    Db for the input grid
 ** \param[in]  dbgrid2    Db for the output grid
 **
 ** \param[out] ind0       Starting address
 ** \param[out] nxyz       Number of subdivisions
 ** \param[out] ntot       Total number of grid nodes
 **
 *****************************************************************************/
static int st_is_subgrid(int verbose,
                         const char *title,
                         DbGrid *dbgrid1,
                         DbGrid *dbgrid2,
                         int *ind0,
                         int *nxyz,
                         int *ntot)
{
  double d;
  int ndim;

  /* Initializations */

  ndim = dbgrid1->getNDim();

  (*ntot) = 1;
  for (int idim = 0; idim < ndim; idim++)
  {
    /* Initializations */

    ind0[idim] = 0;
    nxyz[idim] = 1;

    /* Is origin of the output grid is located on a node of input grid */

    d = (dbgrid2->getX0(idim) - dbgrid1->getX0(idim)) / dbgrid1->getDX(idim);
    if (!isInteger(d))
    {
      messerr(
          "The origin of the Output Grid does not coincide with a node of the Input Grid");
      return (0);
    }
    ind0[idim] = (int) floor(d + 0.5);

    /* Are grid meshes multiple */

    d = dbgrid2->getDX(idim) / dbgrid1->getDX(idim);
    if (!isInteger(d))
    {
      messerr(
          "The grid cell of the Output Grid is not a multiple of the grid cell of the Input Grid");
      return (0);
    }
    nxyz[idim] = (int) floor(d + 0.5);
    (*ntot) *= nxyz[idim];
  }

  if (verbose)
  {
    mestitle(1, title);
    message("- Number of Cells =");
    for (int idim = 0; idim < ndim; idim++)
      message(" %d", nxyz[idim]);
    message("\n");
    message("- Index of Origin =");
    for (int idim = 0; idim < ndim; idim++)
      message(" %d", ind0[idim]);
    message("\n");
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Upscale one variable from a grid Db into another grid Db
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid1    Db for the input grid
 ** \param[in]  dbgrid2    Db for the output grid
 ** \param[in]  orient     Upscaling direction (0 to 2)
 ** \param[in]  verbose    Verbose flag
 **
 *****************************************************************************/
int db_upscale(DbGrid *dbgrid1, DbGrid *dbgrid2, int orient, int verbose)
{
  double *valtab0, *valtab1, *valtab2, *numtab0, *numtab1, *numtab2;
  double result1, result2, result, probtot;
  int error, ndim, ind0[3], nxyz[3], ixyz[3], iech, iptr, ntot, ncol;
  int flag_save, iech_save;

  /* Initializations */

  valtab0 = valtab1 = valtab2 = numtab0 = numtab1 = numtab2 = nullptr;
  error = 1;
  iech_save = (int) get_keypone("Upscale.Converge.Block", 0);

  /* Preliminary checks */

  ndim = dbgrid1->getNDim();
  if (ndim < 1 || ndim > 3)
  {
    messerr("This function is limited to 2-D or 3-D input grids");
    goto label_end;
  }
  if (!dbgrid1->isVariableNumberComparedTo(1)) goto label_end;
  if (orient < 1 || orient > ndim)
  {
    messerr("Inconsistency between Orientation (%d) and Space dimension (%d)",
            orient, ndim);
    goto label_end;
  }
  orient--;

  /* Check that the output grid is a subgrid of the input grid */

  if (!st_is_subgrid(verbose, "Upscaling", dbgrid1, dbgrid2, ind0, nxyz, &ntot))
    goto label_end;

  /* Create the new variable in the output file */

  ncol = 3;
  iptr = dbgrid2->addColumnsByConstant(ncol, TEST);
  if (iptr < 0) goto label_end;
  dbgrid2->setLocatorsByUID(ncol, iptr, ELoc::Z);

  /* Core allocation */

  numtab0 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (numtab0 == nullptr) goto label_end;
  numtab1 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (numtab1 == nullptr) goto label_end;
  numtab2 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (numtab2 == nullptr) goto label_end;
  valtab0 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (valtab0 == nullptr) goto label_end;
  valtab1 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (valtab1 == nullptr) goto label_end;
  valtab2 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (valtab2 == nullptr) goto label_end;

  /* Loop on the cells of the Output Grid */

  for (iech = 0; iech < dbgrid2->getSampleNumber(); iech++)
  {
    result1 = result2 = result = TEST;
    OptDbg::setCurrentIndex(iech + 1);
    flag_save = (iech == iech_save - 1);
    if (dbgrid2->isActive(iech))
    {
      db_index_sample_to_grid(dbgrid2, iech, ixyz);

      /* Load the subgrid to be upscaled */

      probtot = st_extract_subgrid(verbose, 0, iech, dbgrid2->getSampleNumber(),
                                   ntot, dbgrid1, ind0, ixyz, nxyz, numtab0,
                                   valtab0);

      if (probtot > 0)
      {

        /* Upscale the corresponding subgrid of the Input Grid */

        st_upscale(orient, nxyz, flag_save, numtab0, numtab1, numtab2, valtab0,
                   valtab1, valtab2, &result1, &result2);
        result = sqrt(result1 * result2);
      }
      else
      {
        result = result1 = result2 = TEST;
      }
    }

    /* Store the result */

    dbgrid2->setLocVariable(ELoc::Z,iech, 0, result1);
    dbgrid2->setLocVariable(ELoc::Z,iech, 1, result2);
    dbgrid2->setLocVariable(ELoc::Z,iech, 2, result);
  }

  /* Set the error return code */

  error = 0;

  label_end: OptDbg::setCurrentIndex(0);
  mem_free((char* ) numtab0);
  mem_free((char* ) numtab1);
  mem_free((char* ) numtab2);
  mem_free((char* ) valtab0);
  mem_free((char* ) valtab1);
  mem_free((char* ) valtab2);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the euclidean distance between current and initial seed locations
 **
 ** \return Calculated squared distance
 **
 ** \param[in]  orient    Diffusion orientation (0 or the space rank dimension)
 ** \param[in]  ndim      Space dimension
 ** \param[in]  locini    Initial seed positions
 ** \param[in]  loccur    Current seed positions
 **
 *****************************************************************************/
static double st_squared_distance(int orient,
                                  int ndim,
                                  const int *locini,
                                  const int *loccur)
{
  double delta, dist;

  dist = 0.;
  for (int idim = 0; idim < ndim; idim++)
  {
    if (orient != 0 && orient != (idim + 1)) continue;
    delta = locini[idim] - loccur[idim];
    dist += delta * delta;
  }
  return (dist);
}

/****************************************************************************/
/*!
 **  Converts from sample index into grid indices
 **
 ** \param[in]  ndim  Space dimension
 ** \param[in]  ntot  Total number of cells in subgrid
 ** \param[in]  nxyz  Dimensions of the subgrid
 ** \param[in]  iech  Rank of the sample
 **
 ** \param[out]  indg Grid indices
 **
 *****************************************************************************/
static void st_sample_to_grid(int ndim,
                              int ntot,
                              const int *nxyz,
                              int iech,
                              int *indg)
{
  for (int idim = ndim - 1; idim >= 0; idim--)
  {
    ntot /= nxyz[idim];
    indg[idim] = iech / ntot;
    iech -= indg[idim] * ntot;
  }
}

/****************************************************************************/
/*!
 **  Converts from grid indices into sample index
 **
 ** \return The absolute index
 **
 ** \param[in]  ndim  Space dimension
 ** \param[in]  nxyz  Dimensions of the subgrid
 ** \param[in]  indg  Grid indices
 **
 *****************************************************************************/
static int st_grid_to_sample(int ndim, const int *nxyz, const int *indg)
{
  int idim, ival;

  ival = indg[ndim - 1];
  if (ival < 0 || ival >= nxyz[ndim - 1]) return (-1);
  for (idim = ndim - 2; idim >= 0; idim--)
  {
    if (indg[idim] < 0 || indg[idim] >= nxyz[idim]) return (-1);
    ival = ival * nxyz[idim] + indg[idim];
  }
  return (ival);
}

/****************************************************************************/
/*!
 **  Find the location of the cell (within the possible cells) which is the
 **  closest to the target cell provided as argument
 **
 ** \return Rank of the cell
 **
 ** \param[in]  ntot  Number of possibilities
 ** \param[in]  tab   Array containing the sample value
 ** \param[in]  cell  Cell location
 **
 *****************************************************************************/
static int st_fixed_position(int ntot, const double *tab, int cell)
{
  int j;

  for (int i = 0; i < ntot; i++)
  {
    j = cell + i;
    if (j < ntot && tab[j] > 0) return (j);
    j = cell - i;
    if (j >= 0 && tab[j] > 0) return (j);
  }
  return (ntot - 1);
}

/****************************************************************************/
/*!
 **  Find the cell linked to the probability
 **
 ** \return Rank of the cell
 **
 ** \param[in]  ntot  Number of possibilities
 ** \param[in]  tab   Array containing the sample value
 ** \param[in]  proba Local probability
 **
 *****************************************************************************/
static int st_find_cell(int ntot, const double *tab, double proba)
{
  double sum1, sum2;

  sum1 = sum2 = 0.;
  for (int i = 0; i < ntot; i++)
  {
    sum2 += tab[i];
    if (proba >= sum1 && proba < sum2) return (i);
    sum1 = sum2;
  }
  return (ntot - 1);
}

/****************************************************************************/
/*!
 **  Migrate the seed to one of the available neighboring cells
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in]  n_nbgh    Number of neighboring cells
 ** \param[in]  nxyz      Dimensions of the subgrid
 ** \param[in]  nbgh      Array giving the neighboring cell location
 ** \param[in]  valwrk    Working array (Dimension: n_nbgh)
 ** \param[in]  valtab0   Array containing the sample value
 ** \param[in]  locwrk    Working array for shifted seed positions
 ** \param[in,out] loccur    Current seed position
 **
 *****************************************************************************/
static void st_migrate_seed(int ndim,
                            int n_nbgh,
                            int *nxyz,
                            const int *nbgh,
                            double *valwrk,
                            const double *valtab0,
                            int *locwrk,
                            int *loccur)
{
  int iabs, ivois;
  double probtot, proba;

  /* Count the number of available neighboring cells */

  probtot = 0.;
  for (ivois = 0; ivois < n_nbgh; ivois++)
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      locwrk[idim] = loccur[idim] + NBGH(ivois, idim);
      locwrk[idim] = Grid::generateMirrorIndex(nxyz[idim], locwrk[idim]);
    }
    iabs = st_grid_to_sample(ndim, nxyz, locwrk);
    valwrk[ivois] = valtab0[iabs];
    probtot += valwrk[ivois];
  }

  /* Draw a migration at random */

  if (probtot > 0)
  {
    proba = law_uniform(0., probtot);
    ivois = st_find_cell(n_nbgh, valwrk, proba);
    for (int idim = 0; idim < ndim; idim++)
      loccur[idim] += NBGH(ivois, idim);
  }
}

/****************************************************************************/
/*!
 **  Print the position (debug option)
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in]  iseed     Rank of the trajectory (from 0)
 ** \param[in]  iter      Rank of the iteration (from 0)
 ** \param[in]  tab       Array containing the coordinates of the position
 **
 ** \remark  This very verbose option only functions if verbose is TRUE
 ** \remark  and the internal flag DEBUG is TRUE (this requires compiling)
 **
 *****************************************************************************/
static void st_print_position(int ndim, int iseed, int iter, int *tab)
{
  if (!DEBUG) return;
  message("Trajectory %d - Iteration %d:", iseed + 1, iter + 1);
  for (int idim = 0; idim < ndim; idim++)
    message(" %4d", tab[idim]);
  message("\n");
}

/****************************************************************************/
/*!
 **  Calculate the average squared distance as a function of the iteration
 **
 ** \param[in]  orient    Diffusion orientation (0 or the space rank dimension)
 ** \param[in]  ndim      Space dimension
 ** \param[in]  ntot      Total number of cells in the subgrid
 ** \param[in]  nseed     Number of seeds
 ** \param[in]  niter     Number of iterations
 ** \param[in]  n_nbgh    Number of neighboring cells
 ** \param[in]  flag_save If the array must be saved for keypair
 ** \param[in]  probtot   Total probability of joins
 ** \param[in]  nxyz      Dimensions of the subgrid
 ** \param[in]  nbgh      Array giving the neighboring cell location
 ** \param[in]  tabini    Array of initial seed positions
 **                       Dimension: nseed * ndim
 ** \param[in]  tabcur    Array of current seed positions
 **                       Dimension: nseed * ndim
 ** \param[in]  tabwrk    Array of shifted seed positions
 **                       Dimension: ndim
 ** \param[in]  valwrk    Working array
 **                       Dimension: n_nbgh
 ** \param[in]  valtab0   Array containing the sample value
 ** \param[in]  verbose   Verbose option
 **
 ** \param[out] cvdist2   Array of squared distances by iteration
 ** \param[out] trsave    Array of trajectories (saved only if defined)
 **
 ** \remarks If the starting position is fixed, it is specified using:
 **              set_keypair("Fixed_Position",...)
 **
 *****************************************************************************/
static void st_updiff(int orient,
                      int ndim,
                      int ntot,
                      int nseed,
                      int niter,
                      int n_nbgh,
                      int flag_save,
                      double probtot,
                      int *nxyz,
                      int *nbgh,
                      int *tabini,
                      int *tabcur,
                      int *tabwrk,
                      double *valwrk,
                      double *valtab0,
                      int verbose,
                      double *cvdist2,
                      double *trsave)
{
  double d2, dmoy, proba;
  int rank, fixed_position, flag_fixed;

  /* Check if a fixed starting position has been defined */

  fixed_position = (int) get_keypone("Fixed_Position", -1);
  flag_fixed = fixed_position >= 0;

  /* Draw initial seed locations */

  for (int iseed = 0; iseed < nseed; iseed++)
  {
    if (flag_fixed)
    {
      rank = st_fixed_position(ntot, valtab0, fixed_position);
    }
    else
    {
      proba = law_uniform(0., probtot);
      rank = st_find_cell(ntot, valtab0, proba);
    }
    st_sample_to_grid(ndim, ntot, nxyz, rank, &TABINI(iseed, 0));
    for (int idim = 0; idim < ndim; idim++)
      TABCUR(iseed,idim) = TABINI(iseed, idim);
    if (verbose) st_print_position(ndim, iseed, -1, &TABCUR(iseed, 0));
  }

  /* Loop on the iterations */

  for (int iter = 0; iter < niter; iter++)
  {
    dmoy = 0.;

    /* Loop on the seed points */

    for (int iseed = 0; iseed < nseed; iseed++)
    {

      /* Migrate the seed */

      st_migrate_seed(ndim, n_nbgh, nxyz, nbgh, valwrk, valtab0, tabwrk,
                      &TABCUR(iseed, 0));

      /* Optional printout */

      if (verbose) st_print_position(ndim, iseed, iter, &TABCUR(iseed, 0));

      /* Calculate the distance between current and initial positions */

      d2 = st_squared_distance(orient, ndim, &TABINI(iseed, 0),
                               &TABCUR(iseed, 0));

      /* Save the trajectory (optional) */

      if (flag_save && trsave != nullptr)
      {
        for (int idim = 0; idim < ndim; idim++)
          TRAJEC(iseed,iter,idim) = TABCUR(iseed, idim);
      }

      /* Update the mean distance */

      dmoy += d2;
    }
    cvdist2[iter] = dmoy / (double) nseed;
  }
}

/****************************************************************************/
/*!
 **  Update the quantities needed for calculating the linear regression
 **  amongst a set of 2-D points
 **
 ** \param[in]      x         X value
 ** \param[in]      y         Y value
 **
 ** \param[in,out]  count     Number of samples
 ** \param[in,out]  sum_x     Sum of the X variable
 ** \param[in,out]  sum_y     Sum of the Y variable
 ** \param[in,out]  sum_xx    Sum of the X*X variable
 ** \param[in,out]  sum_xy    Sum of the X*Y variable
 **
 *****************************************************************************/
static void st_update_regression(double x,
                                 double y,
                                 double *count,
                                 double *sum_x,
                                 double *sum_y,
                                 double *sum_xx,
                                 double *sum_xy)
{
  (*count) += 1.;
  (*sum_x) += x;
  (*sum_y) += y;
  (*sum_xx) += x * x;
  (*sum_xy) += x * y;
}

/****************************************************************************/
/*!
 **  Derive the diffusion factor from the serie of squared distance as
 **  a function of time
 **
 ** \return The diffusion coefficient
 **
 ** \param[in]  niter      Number of iterations
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  pmid       Percentage of niter to store convergence
 ** \param[in]  flag_save  If the array must be saved for keypair
 ** \param[in]  cvdist2    Array containing the squared distances
 **
 ** \param[out] cvsave     Array containing the storage (Dimension: 3 * niter)
 **
 *****************************************************************************/
static double st_get_diff_coeff(int niter,
                                int verbose,
                                double pmid,
                                int flag_save,
                                double *cvdist2,
                                double *cvsave)
{
  double slope, origin, slope_ref, origin_ref, sum_x, sum_y, sum_xy, sum_xx,
      count;
  double mx, my, var, cov;
  int iter, rank_mid;

  slope_ref = origin_ref = TEST;
  rank_mid = (int) (pmid * niter / 100.);
  count = sum_x = sum_y = sum_xx = sum_xy = origin = slope = 0.;

  for (int jter = 0; jter < niter; jter++)
  {
    iter = niter - jter - 1;

    /* Calculate the average slope */

    st_update_regression((double) (iter + 1), cvdist2[iter], &count, &sum_x,
                         &sum_y, &sum_xx, &sum_xy);
    if (count > 1)
    {
      mx = sum_x / count;
      my = sum_y / count;
      var = sum_xx / count - mx * mx;
      cov = sum_xy / count - mx * my;
      slope = cov / var;
      origin = my - slope * mx;
      if (iter == rank_mid)
      {
        slope_ref = slope;
        origin_ref = origin;
      }
    }

    /* Optional printout */

    if (verbose && !FFFF(slope) && !FFFF(origin))
    {
      message("  Rank=%5d Slope=%lf Origin=%lf (Count=%d)", iter + 1, slope,
              origin, (int) count);
      if (iter == rank_mid) message(" - Stored");
      message("\n");
    }

    /* Optional storage */

    if (flag_save)
    {
      cvsave[3 * iter + 0] = cvdist2[iter];
      cvsave[3 * iter + 1] = slope;
      cvsave[3 * iter + 2] = origin;
    }
  }

  /* Optional printout */

  if (verbose) message("- Slope=%lf Origin=%lf\n", slope_ref, origin_ref);

  /* Optional saving */

  if (flag_save) set_keypair("Diffusion.Converge", 1, niter, 3, cvsave);

  return (slope_ref);
}

/****************************************************************************/
/*!
 **  Calculate the diffusion factor from a grid Db into another grid Db
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid1    Db for the input grid
 ** \param[in]  dbgrid2    Db for the output grid
 ** \param[in]  orient     Diffusion orientation (0 or the space rank dimension)
 ** \param[in]  niter      Number of iterations
 ** \param[in]  nseed      Number of seeds
 ** \param[in]  seed       Seed for the random number generation
 ** \param[in]  verbose    Verbose Option
 **
 ** \remarks The user may specify the type of storage
 ** \remarks 0: average squared distance; 1: slope; 2: origin
 ** \remarks      get.keypair("Diffusion.Converge.Option",0)
 ** \remarks The neighborhood is Block(0) or Cross(1) according to
 ** \remarks      get.keypair("Diffusion.Converge.Morpho",0)
 **
 ** \remarks The rank of the convergence if measured from a Starting to an Ending
 ** \remarks percentages of the 'niter' iterations. The 'Mid' percentage gives
 ** \remarks the value assigned to the output block
 ** \remarks      get.keypair("Diffusion.Converge.PMid",70)
 **
 ** \remarks The user can specify the block of interest by using:
 ** \remarks      get.keypair("Diffusion.Converge.Block")
 ** \remarks For the target block, the array of squared distances as a function
 ** \remarks of the iteration is returned using the keypair mechanism:
 ** \remarks      set.keypair("Diffusion.Converge")
 ** \remarks For the target block, the trajectories are saved using keypair;
 ** \remarks      set.keypair("Diffusion.Trajectory.XX")
 **
 *****************************************************************************/
int db_diffusion(DbGrid *dbgrid1,
                 DbGrid *dbgrid2,
                 int orient,
                 int niter,
                 int nseed,
                 int seed,
                 int verbose)
{
  double *valtab0, *numtab0, *valwrk, *cvdist2, *cvsave, *trsave;
  double diff_coeff, pmid, probtot;
  int error, ndim, ind0[3], nxyz[3], ixyz[3], iech, nech, iptr, opt_center;
  int ntot, iech_save, flag_save, opt_morpho, flag_traj;
  int *tabini, *tabcur, *tabwrk, *numrank, n_nbgh;
  char name[40];
  VectorInt nbgh;

  /* Initializations */

  valtab0 = valwrk = numtab0 = cvdist2 = cvsave = trsave = nullptr;
  tabini = tabcur = tabwrk = numrank = nullptr;
  error = 1;
  iech_save = (int) get_keypone("Diffusion.Converge.Block", 0);
  opt_morpho = (int) get_keypone("Diffusion.Converge.Morpho", 1);
  opt_center = (int) get_keypone("Diffusion.Converge.Center", 1);
  flag_traj = (int) get_keypone("Diffusion.Flag.Trajectory", 0);
  pmid = get_keypone("Diffusion.Converge.PMid", 70.);
  if (seed != 0) law_set_random_seed(seed);

  /* Preliminary checks */

  ndim = dbgrid1->getNDim();
  nech = dbgrid2->getSampleNumber();
  if (ndim < 1 || ndim > 3)
  {
    messerr("This function is limited to 2-D or 3-D input grids");
    goto label_end;
  }
  if (orient != 0 && (orient < 1 || orient > ndim))
  {
    messerr("Argument 'orient' (%d) can be 0 or one of the space dimension",
            orient);
    goto label_end;
  }
  if (!dbgrid1->isVariableNumberComparedTo(1)) goto label_end;
  if (pmid < 5 || pmid > 95)
  {
    messerr("'PMid' must lie between 5% and 95%");
    goto label_end;
  }

  /* Check that the output grid is a subgrid of the input grid */

  if (!st_is_subgrid(verbose, "Diffusion Coefficient", dbgrid1, dbgrid2, ind0,
                     nxyz, &ntot)) goto label_end;

  /* Core allocation */

  tabini = (int*) mem_alloc(sizeof(int) * ndim * nseed, 0);
  if (tabini == nullptr) goto label_end;
  tabcur = (int*) mem_alloc(sizeof(int) * ndim * nseed, 0);
  if (tabcur == nullptr) goto label_end;
  tabwrk = (int*) mem_alloc(sizeof(int) * ndim, 0);
  if (tabwrk == nullptr) goto label_end;
  numrank = (int*) mem_alloc(sizeof(int) * ndim * nseed, 0);
  if (numrank == nullptr) goto label_end;
  numtab0 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (numtab0 == nullptr) goto label_end;
  valtab0 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (valtab0 == nullptr) goto label_end;
  cvdist2 = (double*) mem_alloc(sizeof(double) * niter, 0);
  if (cvdist2 == nullptr) goto label_end;
  cvsave = (double*) mem_alloc(sizeof(double) * niter * 3, 0);
  if (cvsave == nullptr) goto label_end;
  if (flag_traj)
  {
    trsave = (double*) mem_alloc(sizeof(double) * niter * nseed * ndim, 0);
    if (trsave == nullptr) goto label_end;
  }

  /* Allocate the neighboring displacement array */

  nbgh = gridcell_neigh(ndim, opt_morpho, 1, opt_center, verbose);
  n_nbgh = (int) nbgh.size() / ndim;
  valwrk = (double*) mem_alloc(sizeof(double) * n_nbgh, 0);
  if (valwrk == nullptr) goto label_end;

  /* Create the new variable in the output file */

  iptr = dbgrid2->addColumnsByConstant(1, TEST);
  if (iptr < 0) goto label_end;

  /* Loop on the cells of the Output Grid */

  for (iech = 0; iech < nech; iech++)
  {
    diff_coeff = TEST;
    OptDbg::setCurrentIndex(iech + 1);
    flag_save = (iech == iech_save - 1);
    if (dbgrid2->isActive(iech))
    {
      db_index_sample_to_grid(dbgrid2, iech, ixyz);

      /* Load the subgrid to be upscaled */

      probtot = st_extract_subgrid(verbose, 1, iech, nech, ntot, dbgrid1, ind0,
                                   ixyz, nxyz, numtab0, valtab0);

      if (probtot > 0)
      {

        /* Upscale the diffusion */

        st_updiff(orient, ndim, ntot, nseed, niter, n_nbgh, flag_save, probtot,
                  nxyz, nbgh.data(), tabini, tabcur, tabwrk, valwrk, valtab0,
                  verbose, cvdist2, trsave);

        /* Derive the diffusion coefficient */

        diff_coeff = st_get_diff_coeff(niter, verbose, pmid, flag_save, cvdist2,
                                       cvsave);

        /* Save the trajectory (optional) */

        if (flag_save && trsave != nullptr)
        {
          for (int iseed = 0; iseed < nseed; iseed++)
          {
            (void) gslSPrintf(name, "Diffusion.Trajectory.%d", iseed + 1);
            for (int iter = 0; iter < niter; iter++)
              for (int idim = 0; idim < ndim; idim++)
                TRAJEC(iseed,iter,idim) = dbgrid2->getCoordinate(iech, idim) +
                TRAJEC(iseed,iter,idim) * dbgrid1->getDX(idim);
            set_keypair(name, 1, niter, ndim, &TRAJEC(iseed, 0, 0));
          }
        }
      }
    }

    /* Store the result */

    dbgrid2->setArray(iech, iptr, diff_coeff);
  }

  /* Set the error return code */

  error = 0;

  label_end: OptDbg::setCurrentIndex(0);
  mem_free((char* ) tabini);
  mem_free((char* ) tabcur);
  mem_free((char* ) tabwrk);
  mem_free((char* ) numrank);
  mem_free((char* ) valwrk);
  mem_free((char* ) numtab0);
  mem_free((char* ) valtab0);
  mem_free((char* ) cvdist2);
  mem_free((char* ) cvsave);
  mem_free((char* ) trsave);
  return (error);
}

/****************************************************************************/
/*!
 **  Create residuals
 **
 ** \return  Error returned code
 **
 ** \param[in]  verbose Verbose flag
 ** \param[in]  nech    Number of samples
 ** \param[in]  tab     Array of sample values (Dimension: nech)
 ** \param[in]  ncut    Number of cutoffs
 ** \param[in]  zcut    Array of cutoff values (Dimension: ncut)
 **
 ** \param[out] nsorted   Number of sorted samples
 ** \param[out] mean      Average of the active data
 ** \param[out] residuals Array of residuals (Dimension: ncut * nech)
 ** \param[out] T         Array of for tonnage
 ** \param[out] Q         Array of for metal quantity
 **
 *****************************************************************************/
int stats_residuals(int verbose,
                    int nech,
                    const double *tab,
                    int ncut,
                    double *zcut,
                    int *nsorted,
                    double *mean,
                    double *residuals,
                    double *T,
                    double *Q)
{
  double value, moyenne;
  int iech, icut, jcut, nactive;

  /* Initializations */

  nactive = (*nsorted) = 0;
  moyenne = 0.;
  for (icut = 0; icut < ncut; icut++)
  {
    T[icut] = Q[icut] = 0.;
    for (iech = 0; iech < nech; iech++)
      RESIDUALS(icut,iech) = 0.;
  }

  /* Loop on the samples to calculate the indicators */

  for (iech = 0; iech < nech; iech++)
  {
    value = tab[iech];
    if (FFFF(value)) continue;
    moyenne += value;
    nactive++;

    /* Loop on the cutoffs */

    for (icut = 0; icut < ncut; icut++)
    {
      if (value < zcut[icut]) continue;
      RESIDUALS(icut,iech) = 1.;
      Q[icut] += value;
      T[icut] += 1.;
    }
  }
  if (nactive <= 0)
  {
    messerr("The calculation failed as there is no active sample");
    return (1);
  }

  /* Calculate the tonnage and meal quantity per class */

  moyenne /= (double) nactive;
  for (icut = 0; icut < ncut; icut++)
  {
    T[icut] /= (double) nactive;
    Q[icut] /= (double) nactive;
  }

  /* Calculate the residuals */

  for (iech = 0; iech < nech; iech++)
  {
    value = tab[iech];
    if (FFFF(value)) continue;

    /* Loop on the cutoffs */

    for (icut = ncut - 1; icut >= 0; icut--)
    {
      value = RESIDUALS(icut,iech) / T[icut];
      if (icut > 0)
      {
        jcut = icut - 1;
        value -= RESIDUALS(jcut,iech) / T[jcut];
      }
      else
      {
        value -= 1.;
      }
      RESIDUALS(icut,iech) = value;
    }
  }

  /* Verbose optional option */

  if (verbose)
  {
    mestitle(0, "Building residuals");
    message("Number of sorted samples = %d\n", nactive);
    for (icut = 0; icut < ncut; icut++)
      message("Cutoff %2d (above %lf) - Tonnage = %lf - Metal = %lf\n",
              icut + 1, zcut[icut], T[icut], Q[icut]);
  }

  (*nsorted) = nactive;
  (*mean) = moyenne;
  return (0);
}
