/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "geoslib_old_f.h"
#include "Basic/Utilities.hpp"
#include "Db/Db.hpp"

#include <math.h>

typedef struct
{
  double xg;
  double yg;
  double zg;
  double extx;
  double exty;
  double vx[3];
  double vy[3];
  double vz[3];
} Surf_Def;

#define MAT(i,j)    (mat[3*(i) + (j)])
#define EIGVEC(i,j) (eigvec[3*(i) + (j)])

static int VERBOSE = 0;

/*****************************************************************************/
/*!
 **  Manage the Projection system
 **
 ** \return  Pointer to the newly managed Surf_Def structure
 **
 ** \param[in]  mode           Type operation (1: allocation; -1: deallocation)
 ** \param[in]  surf_reference Surf_Def structure (used for deallocation)
 **
 *****************************************************************************/
static Surf_Def* st_reference_manage(int mode, Surf_Def *surf_reference)
{
  /* Dispatch */

  if (mode > 0)
  {

    /* Allocation */

    surf_reference = (Surf_Def*) mem_alloc(sizeof(Surf_Def), 1);
    surf_reference->xg = 0.;
    surf_reference->yg = 0.;
    surf_reference->zg = 0.;
    surf_reference->extx = 0.;
    surf_reference->exty = 0.;
    for (int i = 0; i < 3; i++)
    {
      surf_reference->vx[i] = 0.;
      surf_reference->vy[i] = 0.;
      surf_reference->vz[i] = 0.;
    }
  }
  else
  {

    /* Deallocation */

    surf_reference = (Surf_Def*) mem_free((char* ) surf_reference);
  }
  return (surf_reference);
}

/*****************************************************************************/
/*!
 **  Normalize the input vector
 **
 ** \param[in,out]  vect    Vector to be normalized
 **
 *****************************************************************************/
static void st_normalize_vector(double *vect)
{
  double norme;

  norme = 0;
  for (int i = 0; i < 3; i++)
    norme += vect[i] * vect[i];

  if (norme <= 0) return;
  norme = sqrt(norme);

  for (int i = 0; i < 3; i++)
    vect[i] /= norme;
}

/*****************************************************************************/
/*!
 **  Define the Projection new system by Total Least Squares
 **
 ** \param[in]  db             Db structure
 ** \param[in]  iptr_init      Array of coordinate locators (initial)
 **
 ** \param[out] surf_reference Surf_Def structure
 **
 *****************************************************************************/
static int st_reference_define(Db *db, int *iptr_init, Surf_Def *surf_reference)
{
  double mat[9], eigval[3], eigvec[9], gg[3], nn, x, y, z;
  int ix, iy, iz;

  /* Calculate the center of gravity */

  nn = 0;
  for (int i = 0; i < 3; i++)
    gg[i] = 0.;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    x = db->getArray(iech, iptr_init[0]);
    y = db->getArray(iech, iptr_init[1]);
    z = db->getArray(iech, iptr_init[2]);
    gg[0] += x;
    gg[1] += y;
    gg[2] += z;
    nn += 1.;
  }

  surf_reference->xg = gg[0] / nn;
  surf_reference->yg = gg[1] / nn;
  surf_reference->zg = gg[2] / nn;

  /* Establish the mean plane */

  for (int i = 0; i < 9; i++)
    mat[i] = 0.;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;
    x = db->getArray(iech, iptr_init[0]) - surf_reference->xg;
    y = db->getArray(iech, iptr_init[1]) - surf_reference->yg;
    z = db->getArray(iech, iptr_init[2]) - surf_reference->zg;
    MAT(0,0) += x * x;
    MAT(0,1) += x * y;
    MAT(0,2) += x * z;
    MAT(1,0) += y * x;
    MAT(1,1) += y * y;
    MAT(1,2) += y * z;
    MAT(2,0) += z * x;
    MAT(2,1) += z * y;
    MAT(2,2) += z * z;
  }

  /* Eigen values decomposition */

  if (matrix_eigen(mat, 3, eigval, eigvec))
  {
    messerr("Error in the Plane determination");
    return (1);
  }

  /* Look for the smallest eigen value */

  iz = 0;
  for (int i = 1; i < 3; i++)
    if (eigval[i] < eigval[iz]) iz = i;
  ix = (iz >= 2) ? 0 :
                   iz + 1;
  iy = 3 - ix - iz;

  /* Calculate the direction vectors of the reference system */

  surf_reference->vx[0] = EIGVEC(ix, 0);
  surf_reference->vx[1] = EIGVEC(ix, 1);
  surf_reference->vx[2] = EIGVEC(ix, 2);
  st_normalize_vector(surf_reference->vx);

  surf_reference->vy[0] = EIGVEC(iy, 0);
  surf_reference->vy[1] = EIGVEC(iy, 1);
  surf_reference->vy[2] = EIGVEC(iy, 2);
  st_normalize_vector(surf_reference->vy);

  surf_reference->vz[0] = EIGVEC(iz, 0);
  surf_reference->vz[1] = EIGVEC(iz, 1);
  surf_reference->vz[2] = EIGVEC(iz, 2);
  st_normalize_vector(surf_reference->vz);

  if (VERBOSE)
  {
    print_matrix("Gravity Center ", 0, 1, 3, 1, NULL, gg);
    print_matrix("Eigen Values   ", 0, 1, 3, 1, NULL, eigval);
    print_matrix("Eigen Vector #1", 0, 1, 3, 1, NULL, surf_reference->vx);
    print_matrix("Eigen Vector #2", 0, 1, 3, 1, NULL, surf_reference->vy);
    print_matrix("Eigen Vector #3", 0, 1, 3, 1, NULL, surf_reference->vz);
  }

  return (0);
}

/*****************************************************************************/
/*!
 **  Initialize the minimum-maximum array
 **
 ** \param[out]  mima      Statistics array
 ** \param[out]  nundefs   Number of undefined values
 **
 *****************************************************************************/
static void st_mima_init(double mima[2][3], int *nundefs)
{
  (*nundefs) = 0;
  for (int idim = 0; idim < 3; idim++)
  {
    mima[0][idim] = +1.e30;
    mima[1][idim] = -1.e30;
  }
}

/*****************************************************************************/
/*!
 **  Print the minimum-maximum array
 **
 ** \param[in,out] mima    Statistics array
 ** \param[in]     x,y,z   Coordinates
 ** \param[in,out] nundefs Number of undefined values
 **
 *****************************************************************************/
static void st_mima_process(double mima[2][3],
                            double x,
                            double y,
                            double z,
                            int *nundefs)
{
  if (FFFF(x)) (*nundefs)++;
  if (x < mima[0][0]) mima[0][0] = x;
  if (x > mima[1][0]) mima[1][0] = x;
  if (y < mima[0][1]) mima[0][1] = y;
  if (y > mima[1][1]) mima[1][1] = y;
  if (z < mima[0][2]) mima[0][2] = z;
  if (z > mima[1][2]) mima[1][2] = z;
}

/*****************************************************************************/
/*!
 **  Update the minimum-maximum array
 **
 ** \param[in]  title     Title of the statistics
 ** \param[in]  mima      Statistics array
 ** \param[in]  nundefs   Number of undefined values
 **
 *****************************************************************************/
static void st_mima_print(const char *title, double mima[2][3], int nundefs)
{
  message("Statistics for %s\n", title);
  if (nundefs > 0) message("- Number of Undefined values = %d\n", nundefs);
  message("- X : %15.4lf %15.4lf - Delta : %lf\n", mima[0][0], mima[1][0],
          mima[1][0] - mima[0][0]);
  message("- Y : %15.4lf %15.4lf - Delta : %lf\n", mima[0][1], mima[1][1],
          mima[1][1] - mima[0][1]);
  message("- Z : %15.4lf %15.4lf - Delta : %lf\n", mima[0][2], mima[1][2],
          mima[1][2] - mima[0][2]);
}

/*****************************************************************************/
/*!
 **  Project the samples on the reference plane
 **
 ** \param[in]  surf_reference Surf_Ref structure
 ** \param[in]  db             Db structure
 ** \param[in]  iptr_init      Array of coordinate locators (initial)
 ** \param[in]  iptr_proj      Array of coordinate locators (projected)
 **
 *****************************************************************************/
static void st_transform_init2proj(Surf_Def *surf_reference,
                                   Db *db,
                                   int *iptr_init,
                                   int *iptr_proj)
{
  double x, y, z, newx, newy, newz, mima_in[2][3], mima_out[2][3];
  int n, nundefs_in, nundefs_out;

  st_mima_init(mima_in, &nundefs_in);
  st_mima_init(mima_out, &nundefs_out);

  n = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Load the initial coordinates */

    x = db->getArray(iech, iptr_init[0]);
    y = db->getArray(iech, iptr_init[1]);
    z = db->getArray(iech, iptr_init[2]);

    /* Update statistics */

    st_mima_process(mima_in, x, y, z, &nundefs_in);

    /* Translate to the origin */

    x -= surf_reference->xg;
    y -= surf_reference->yg;
    z -= surf_reference->zg;

    /* Calculate the new coordinates */

    newx = (x * surf_reference->vx[0] + y * surf_reference->vx[1]
            + z * surf_reference->vx[2]);
    newy = (x * surf_reference->vy[0] + y * surf_reference->vy[1]
            + z * surf_reference->vy[2]);
    newz = (x * surf_reference->vz[0] + y * surf_reference->vz[1]
            + z * surf_reference->vz[2]);

    /* Update statistics */

    st_mima_process(mima_out, newx, newy, newz, &nundefs_out);

    /* Store the new coordinates */

    db->setArray(iech, iptr_proj[0], newx);
    db->setArray(iech, iptr_proj[1], newy);
    db->setArray(iech, iptr_proj[2], newz);
    n++;
  }

  /* Store the rectangle extension */

  surf_reference->extx = (mima_out[1][0] - mima_out[0][0]) / 2.;
  surf_reference->exty = (mima_out[1][1] - mima_out[0][1]) / 2.;

  if (VERBOSE && n > 0)
  {
    message("From Initial to Projected Systems (%d values)\n", n);
    st_mima_print("Initial System", mima_in, nundefs_in);
    st_mima_print("Projected System", mima_out, nundefs_out);
  }
}

/*****************************************************************************/
/*!
 **  Project the samples from the projected to the initial space
 **
 ** \param[in]  surf_reference Surf_Ref structure
 ** \param[in]  npoint         Number of points
 ** \param[in]  points         Array of 3-D coordinates
 **
 *****************************************************************************/
static void st_transform_proj2init(Surf_Def *surf_reference,
                                   int npoint,
                                   VectorDouble& points)
{
  double x, y, z, newx, newy, newz, mima_in[2][3], mima_out[2][3];
  int n, nundefs_in, nundefs_out;

  st_mima_init(mima_in, &nundefs_in);
  st_mima_init(mima_out, &nundefs_out);

  n = 0;
  for (int ip = 0; ip < npoint; ip++)
  {

    /* Load the initial coordinates */

    x = points[3 * ip];
    y = points[3 * ip + 1];
    z = points[3 * ip + 2];

    /* Update statistics */

    st_mima_process(mima_in, x, y, z, &nundefs_in);

    /* Calculate the new coordinates */

    newx = (x * surf_reference->vx[0] + y * surf_reference->vy[0]
            + z * surf_reference->vz[0]);
    newy = (x * surf_reference->vx[1] + y * surf_reference->vy[1]
            + z * surf_reference->vz[1]);
    newz = (x * surf_reference->vx[2] + y * surf_reference->vy[2]
            + z * surf_reference->vz[2]);

    newx += surf_reference->xg;
    newy += surf_reference->yg;
    newz += surf_reference->zg;

    /* Update statistics */

    st_mima_process(mima_out, newx, newy, newz, &nundefs_out);

    /* Store the new coordinates */

    points[3 * ip] = newx;
    points[3 * ip + 1] = newy;
    points[3 * ip + 2] = newz;
    n++;
  }

  if (VERBOSE && n > 0)
  {
    message("From Projected to Initial Systems (%d values)\n", n);
    st_mima_print("Projected System", mima_in, nundefs_in);
    st_mima_print("Initial System", mima_out, nundefs_out);
  }
}

/*****************************************************************************/
/*!
 **  Select the samples corresponding to the target code (if defined)
 **
 ** \return  Returns the number of selected samples
 **
 ** \param[in]  db             Db structure
 ** \param[in]  icode          Reference Code value
 ** \param[in]  iptr_sel       Pointer to the selection
 **
 ** \remarks Any already existing selection is not taken into account
 ** \remarks If the code is undefined, this function has no action
 **
 *****************************************************************************/
static int st_selection_per_code(Db *db, int icode, int iptr_sel)
{
  int number;

  number = 0;
  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (IFFFF(icode) || (int) db->getLocVariable(ELoc::C,iech,0) == icode)
    {

      db->setArray(iech, iptr_sel, 1.);
      number++;
    }
    else
      db->setArray(iech, iptr_sel, 0.);
  }
  return (number);
}

/*****************************************************************************/
/*!
 **  Concatenate the resulting arrays to the returned arrays
 **
 ** \param[in]  ndim         Space dimension
 ** \param[in]  ntri         Number of triangles
 ** \param[in]  npoints      Number of vertices
 ** \param[in]  triloc       Array on the triangle corners
 ** \param[in]  poiloc       Array on the vertices coordinates
 **
 ** \param[out] ntri_arg     Cumulated number of triangles
 ** \param[out] npoint_arg   Cumulated number of vertices
 ** \param[out] triangles    Cumulated array on the triangle corners
 ** \param[out] points       Cumulated array on the 3-D vertices coordinates
 **
 *****************************************************************************/
static int st_concatenate_arrays(int ndim,
                                 int ntri,
                                 int npoints,
                                 VectorInt& triloc,
                                 VectorDouble& poiloc,
                                 int *ntri_arg,
                                 int *npoint_arg,
                                 VectorInt& triangles,
                                 VectorDouble& points)
{
  int ecr;

  /* Initial assignments */

  int ntloc = *ntri_arg;
  int nploc = *npoint_arg;

  /* Concatenate the triangles */

  if (ntri > 0)
  {
    int newnum = ntloc + ntri;
    triangles.resize(3 * newnum, 0);

    ecr = ntloc * 3;
    for (int i = 0; i < ntri * 3; i++)
      triangles[ecr++] = triloc[i] + nploc;

    ntloc = newnum;
  }

  /* Concatenate the points */

  if (npoints > 0)
  {
    int newnum = nploc + npoints;
    points.resize(ndim * newnum, 0);

    ecr = nploc * ndim;
    for (int i = 0; i < npoints * ndim; i++)
      points[ecr++] = poiloc[i];

    nploc = newnum;
  }

  /* Set the returned arguments */

  *ntri_arg = ntloc;
  *npoint_arg = nploc;

  return (0);
}

/*****************************************************************************/
/*!
 **  Generate the surface as the rectangle containing the surface
 **
 ** \param[in]  surf_reference Surf_Def structure (used for deallocation)
 **
 ** \param[out] ntri_arg     Number of triangles generated
 ** \param[out] npoint_arg   Number of vertices
 ** \param[out] triangles    Array of triangulate vertex indices
 ** \param[out] points       Array containing the 2-D vertices
 **
 *****************************************************************************/
static int st_rectangle_surface(Surf_Def *surf_reference,
                                int *ntri_arg,
                                int *npoint_arg,
                                VectorInt& triangles,
                                VectorDouble& points)
{
  int ecr;
  static double ratio = 1.5;

  /* Initializations */

  int ndim = 3;
  int ncoord = 5;
  int ntri = 4;
  *ntri_arg = *npoint_arg = 0;

  /* Core allocation */

  triangles.resize(3 * ntri, 0);
  points.resize(ndim * ncoord, 0);

  /* Load the points */

  ecr = 0;
  points[ecr++] = 0.;
  points[ecr++] = 0.;
  points[ecr++] = 0.;

  points[ecr++] = ratio * surf_reference->extx;
  points[ecr++] = ratio * surf_reference->exty;
  points[ecr++] = 0.;

  points[ecr++] = -ratio * surf_reference->extx;
  points[ecr++] = ratio * surf_reference->exty;
  points[ecr++] = 0.;

  points[ecr++] = -ratio * surf_reference->extx;
  points[ecr++] = -ratio * surf_reference->exty;
  points[ecr++] = 0.;

  points[ecr++] = ratio * surf_reference->extx;
  points[ecr++] = -ratio * surf_reference->exty;
  points[ecr++] = 0.;

  /* Load the triangles */

  ecr = 0;
  triangles[ecr++] = 1;
  triangles[ecr++] = 2;
  triangles[ecr++] = 3;

  triangles[ecr++] = 1;
  triangles[ecr++] = 3;
  triangles[ecr++] = 4;

  triangles[ecr++] = 1;
  triangles[ecr++] = 4;
  triangles[ecr++] = 5;

  triangles[ecr++] = 1;
  triangles[ecr++] = 5;
  triangles[ecr++] = 2;

  /* Set the error code */

  *ntri_arg = ntri;
  *npoint_arg = ncoord;
  return 0;
}

/*****************************************************************************/
/*!
 **  Free the triangleio structure
 **
 ** \param[in]  db        Db structure
 ** \param[in]  model     Model structure
 ** \param[in]  triswitch Triangulation option
 ** \param[in]  icode0    Reference Code attributed to the Target Fault
 ** \param[in]  verbose   Verbose option
 **
 ** \param[out] ncode_arg     Number of different codes
 ** \param[out] ntri_arg      Number of triangles
 ** \param[out] npoint_arg    Number of vertices
 ** \param[out] codesel       Selected code (if any)
 ** \param[out] ntcode        Array for the number of triangles per code
 ** \param[out] triangles     Array on the triangle corners
 ** \param[out] points        Array on the 3-D vertices coordinates
 **
 ** \remarks The returned arrays 'triangle', 'points'
 ** \remarks must be freed by the calling function
 **
 *****************************************************************************/
int db_trisurf(Db *db,
               Model *model,
               const String &triswitch,
               int icode0,
               int verbose,
               int *ncode_arg,
               int *ntri_arg,
               int *npoint_arg,
               double *codesel,
               VectorInt& ntcode,
               VectorInt& triangles,
               VectorDouble& points)
{
  Surf_Def *surf_reference;
  int iptr_sel, iptr_init[3], iptr_proj[3], error;
  int ndim, ntriloc, ntricum, npoicum, npoiloc, icode, ncodes, ncode_eff, number;
  int flag_rectangle_surface;
  VectorDouble codetab;
  SPDE_Option s_option;
  VectorInt triloc;
  VectorDouble poiloc;

  /* Initializations */

  error = 1;
  ntricum = npoicum = *ncode_arg = *ntri_arg = *npoint_arg = 0;
  VERBOSE = verbose;
  ndim = db->getNDim();
  iptr_sel = -1;
  surf_reference = (Surf_Def*) NULL;
  for (int idim = 0; idim < ndim; idim++)
    iptr_init[idim] = iptr_proj[idim] = -1;
  flag_rectangle_surface = (int) get_keypone("Flag_Rectangle_Surface", 0);

  /* Process the code */

  ncode_eff = 1;
  ncodes = 1;
  if (db->hasLocVariable(ELoc::C))
  {
    codetab = db->getCodeList();
    ncode_eff = static_cast<int>(codetab.size());
    if (!IFFFF(icode0))
    {
      if (icode0 < 0 || icode0 >= ncode_eff)
      {
        messerr("'code'(%d) is incompatible with total number of codes (%d)",
                icode0, ncode_eff);
        goto label_end;
      }
    }
    else
    {
      ncodes = ncode_eff;
    }
  }

  /* Preliminary checks */

  if (ndim != 3)
  {
    messerr("This function is restricted to the case of 3-D Db");
    goto label_end;
  }
  s_option = spde_option_alloc();
  spde_option_update(s_option, triswitch);

  /* Allocate the reference projection */

  surf_reference = st_reference_manage(1, NULL);
  if (surf_reference == (Surf_Def*) NULL) goto label_end;

  /* Core allocation */

  ntcode.resize(ncode_eff, 0);

  /* Create the new variables */

  for (int idim = 0; idim < ndim; idim++)
  {
    iptr_init[idim] = db->getColIdxByLocator(ELoc::X, idim);
    iptr_proj[idim] = db->addColumnsByConstant(1, TEST);
    if (iptr_proj[idim] < 0) goto label_end;
  }
  iptr_sel = db->addColumnsByConstant(1, 1.);
  if (iptr_sel < 0) goto label_end;
  db->setLocatorByUID(iptr_sel, ELoc::SEL);

  /* Set the new 2-D coordinates and turn the third coordinate into variable */

  db->clearLocators(ELoc::X);
  for (int idim = 0; idim < 2; idim++)
    db->setLocatorByUID(iptr_proj[idim], ELoc::X, idim);
  db->setLocatorByUID(iptr_proj[ndim - 1], ELoc::Z);

  /* Loop on the set of points per Fault */

  for (int jcode = 0; jcode < ncodes; jcode++)
  {
    ntriloc = npoiloc = 0;
    if (db->hasLocVariable(ELoc::C))
    {
      icode = (!IFFFF(icode0)) ? (int) codetab[icode0] : (int) codetab[jcode];
      message("\nProcessing Fault for code %d\n", icode);
    }
    else
      icode = ITEST;

    /* Create the selection (optional) */

    number = st_selection_per_code(db, icode, iptr_sel);
    if (number <= 0) goto label_suite;

    /* Define the reference projection surface */

    if (st_reference_define(db, iptr_init, surf_reference)) goto label_end;

    /* Project the data */

    st_transform_init2proj(surf_reference, db, iptr_init, iptr_proj);

    if (! flag_rectangle_surface)
    {

      /* Perform the estimation of the elevation (in the projected space) */

      if (kriging2D_spde(db, model, s_option, 0, &ntriloc, &npoiloc, triloc,
                         poiloc)) goto label_end;
    }
    else
    {

      /* Generate the surface as the rectangle containing the fault */

      if (st_rectangle_surface(surf_reference, &ntriloc, &npoiloc, triloc,
                               poiloc)) goto label_end;
    }

    /* Unproject the results */

    st_transform_proj2init(surf_reference, npoiloc, poiloc);

    /* Concatenate the resulting arrays */

    if (st_concatenate_arrays(3, ntriloc, npoiloc, triloc, poiloc, &ntricum,
                              &npoicum, triangles, points)) goto label_end;

    /* Verbose output */

    message("- Count of samples   = %d\n", number);
    message("- Count of triangles = %d\n", ntriloc);

    label_suite:
    ntcode[jcode] = ntriloc;
  }

  /* Set the error return code */

  (*ncode_arg) = ncodes;
  *ntri_arg = ntricum;
  *npoint_arg = npoicum;
  *codesel = (db->hasLocVariable(ELoc::C) && !IFFFF(icode0)) ? codetab[icode0] : TEST;
  error = 0;

  label_end:

  /* Free memory */

  surf_reference = st_reference_manage(-1, surf_reference);

  /* Delete new temporary variables */

  db->deleteColumnByUID(iptr_sel);
  for (int idim = 0; idim < ndim; idim++)
    db->deleteColumnByUID(iptr_proj[ndim - idim - 1]);

  return (error);
}
