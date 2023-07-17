/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"

#include "Enum/EMorpho.hpp"

#include "Basic/Utilities.hpp"
#include "Db/Db.hpp"
#include "Morpho/Morpho.hpp"

#include <math.h>

static int RADIUS[3];
static int LARGE = 9999999;


/*! \cond */
#define CROSS 0
#define BLOCK 1

#define IMAGES(ix,iy,iz) ((option == BLOCK) ? imagtmp.getValue(ix,iy,iz) : imagin.getValue(ix,iy,iz))

#define GRID_ADD(ix,iy,iz)  ((ix)+(NX[0]*((iy)+NX[1]*(iz))))
#define TAB(i,j,k)          (tab[GRID_ADD(i,j,k)])
#define COMPNUM(i,j,k)      (compnum[GRID_ADD(i,j,k)])
/*! \endcond */

/*****************************************************************************/
/*!
 **  Defines the image radius (global variables)
 **
 ** \param[in]  radius Radius of the structural element
 **
 *****************************************************************************/
void _st_morpho_image_radius_define(const VectorInt &radius)
{
  int size = static_cast<int>(radius.size());
  RADIUS[0] = (size > 0) ? radius[0] : 0;
  RADIUS[1] = (size > 1) ? radius[1] : 0;
  RADIUS[2] = (size > 2) ? radius[2] : 0;
}

/*****************************************************************************/
/*!
 **  Returns the sizes of all the connected components
 **
 ** \return  Total volume of the measured connected components
 **
 ** \param[in]  compnum array containing the component index
 ** \param[in]  nbcomp  number of connected components to be measured
 **
 ** \param[out] sizes array containing the sizes of the connected components
 **                   labelled from 1 to nbcomp
 **
 *****************************************************************************/
int _st_morpho_label_size(const VectorDouble &compnum,
                          int nbcomp,
                          VectorInt &sizes)
{
  int total = 0;
  for (int i = 0; i < (int) compnum.size(); i++)
  {
    int val = (int) compnum[i];
    if (val > 0 && val <= nbcomp)
    {
      sizes[val - 1]++;
      total++;
    }
  }
  return (total);
}

/*****************************************************************************/
/*!
 **  Orders the connected components using a rank array
 **
 ** \param[in]  compnum array containing the component index
 ** \param[in]  order   array containing the ordering criterion
 ** \param[in]  nbcomp  number of connected components
 **
 *****************************************************************************/
void _st_morpho_label_order(VectorDouble &compnum,
                            const VectorInt &order,
                            int nbcomp)
{
  for (int i = 0; i < (int) compnum.size(); i++)
  {
    int val = (int) compnum[i];
    if (val <= 0) continue;
    int found = -1;
    for (int j = nbcomp - 1; j >= 0 && found < 0; j--)
      if (val == order[j]) found = j;
    if (found < 0) messageAbort("st_morpho_label_order");
    compnum[i] = nbcomp - found;
  }
  return;
}

/**
 * Duplicates the contenst of the input BImage into the output BImage
 * @param imagin
 * @param imagout
 */
void morpho_duplicate(const BImage &imagin, BImage &imagout)
{
  imagout = imagin;
}

/**
 * Labels the connected components for a BImage
 *
 * @remarks The labels are sorted by decreasing sizes and an optional message
 * is issued for displaying the component sizes
 */
VectorDouble morpho_labelling(int option,
                              int flag_size,
                              const BImage &imagin,
                              double ccvoid,
                              bool verbose)
{
  int jx, jy, jz, total, count, ref, iad, local, part_grain;
  int i, size, nbtest, nbcomp, ix, iy, iz, icomp, ival, itest, il, jcomp[26];
  VectorInt list_array, sizes, order;
  VectorDouble compnum;
  int id[26][3] = { { -1, 00, 00 },
                    { 00, -1, 00 },
                    { 00, 00, -1 },
                    { 01, 00, 00 },
                    { 00, 01, 00 },
                    { 00, 00, 01 },
                    { -1, -1, -1 },
                    { 00, -1, -1 },
                    { 01, -1, -1 },
                    { -1, 00, -1 },
                    { -1, 01, -1 },
                    { 00, 01, -1 },
                    { 01, 00, -1 },
                    { 01, 01, -1 },
                    { -1, -1, 00 },
                    { -1, 01, 00 },
                    { 01, -1, 00 },
                    { 01, 01, 00 },
                    { -1, -1, 01 },
                    { 00, -1, 01 },
                    { 01, -1, 01 },
                    { -1, 00, 01 },
                    { -1, 01, 01 },
                    { 00, 01, 01 },
                    { 01, 00, 01 },
                    { 01, 01, 01 } };
  int ndel[2] = { 6, 26 };
  int quantum = 100;

  /* Initializations */

  nbtest = ndel[option];
  nbcomp = total = 0;
  VectorInt NX = imagin.getNDimsExt(3);

  /* Attempt to allocate the initial quantum */

  int nxyz = imagin.getNPixels();
  size = MIN(quantum, nxyz);
  list_array.resize(size);
  compnum.resize(nxyz);
  for (i = 0; i < nxyz; i++) compnum[i] = 0.;

  /* Numbering the connected components */

  for (iz = 0; iz < NX[2]; iz++)
    for (iy = 0; iy < NX[1]; iy++)
      for (ix = 0; ix < NX[0]; ix++)
      {
        if (imagin.getValue(ix, iy, iz))
        {
          icomp = LARGE;
          for (itest = 0; itest < nbtest; itest++)
          {
            jcomp[itest] = LARGE;
            jx = ix + id[itest][0];
            jy = iy + id[itest][1];
            jz = iz + id[itest][2];
            if (imagin.isInside(jx, jy, jz))
            {
              ival = (int) COMPNUM(jx, jy, jz);
              if (ival > 0)
              {
                jcomp[itest] = list_array[ival - 1];
                if (jcomp[itest] < icomp) icomp = jcomp[itest];
              }
            }
          }
          if (icomp != LARGE)
          {
            for (itest = 0; itest < nbtest; itest++)
              if (jcomp[itest] != LARGE && jcomp[itest] != icomp)
                for (il = jcomp[itest]; il <= nbcomp; il++)
                  if (list_array[il - 1] == jcomp[itest])
                    list_array[il - 1] = icomp;
          }
          else
          {
            nbcomp++;
            if (nbcomp > size)
            {
              size += quantum;
              list_array.resize(size);
            }
            list_array[nbcomp - 1] = nbcomp;
            icomp = nbcomp;
          }
        }
        else
        {
          icomp = 0;
        }
        COMPNUM(ix,iy,iz) = (double) icomp;
      }

  /* Compressing the list */

  for (il = ival = 0; il < nbcomp; il++)
    if (list_array[il] == il + 1)
      list_array[il] = ++ival;
    else
      list_array[il] = list_array[list_array[il] - 1];
  nbcomp = ival;

  /* Update the planes */

  for (iz = 0; iz < imagin.getNDims(2); iz++)
    for (iy = 0; iy < imagin.getNDims(1); iy++)
      for (ix = 0; ix < imagin.getNDims(0); ix++)
      {
        iad = (int) COMPNUM(ix, iy, iz);
        if (iad != 0)
          COMPNUM(ix,iy,iz) = list_array[iad - 1];
        else
          COMPNUM(ix,iy,iz) = TEST;
      }

  /* Order the components */

  if (nbcomp > 0)
  {
    sizes.resize(nbcomp, 0);
    total = _st_morpho_label_size(compnum, nbcomp, sizes);
    order.resize(nbcomp);
    for (i = 0; i < nbcomp; i++)
      order[i] = i + 1;
    VH::arrangeInPlace(1, order, sizes, true, nbcomp);
    _st_morpho_label_order(compnum, order, nbcomp);

    if (flag_size)
      for (i = 0; i < nxyz; i++)
      {
        if (FFFF(compnum[i])) continue;
        ival = (int) compnum[i];
        compnum[i] = sizes[ival - 1];
      }

    for (i = 0; i < nxyz; i++)
      if (FFFF(compnum[i])) compnum[i] = ccvoid;
  }

  /* Display the label component size statistics */

  if (verbose)
  {
    if (nbcomp == 0)
    {
      message("No grain has been detected\n");
    }
    else
    {
      message("Number of connected components = %d\n", nbcomp);
      message("Total connected volume = %d\n",total);
      message("   Component     Number         Total     Cumul (percent)\n");
      ref = sizes[order[nbcomp - 1] - 1];
      count = part_grain = 0;
      for (i = nbcomp - 1; i >= 0; i--)
      {
        local = sizes[order[i] - 1];

        if (local != ref)
        {
          message("%12d  %9d  %12d      %7.3f\n", ref, count, count * ref,
                  100. * part_grain / total);
          ref = local;
          count = 0;
        }

        part_grain += local;
        count++;
      }
      message("%12d  %9d  %12d      %7.3f\n", ref, count, count * ref,
              100. * part_grain / total);
    }
  }
  return compnum;
}

/**
 * Returns the array of dimensions of the connex components
 * @remarks The labels are sorted by decreasing sizes and an optional message
 * is issued for displaying the component sizes
 */
VectorInt morpho_labelsize(int option, const BImage& imagin)
{
  VectorInt sizes;
  VectorDouble compnum = morpho_labelling(option, 0, imagin, TEST);
  int nbcomp = compnum.size();
  if (nbcomp > 0)
  {
    sizes.resize(nbcomp, 0);
    (void) _st_morpho_label_size(compnum, nbcomp, sizes);
  }
  return sizes;
}

/**
 * Perform the "erosion" of the input image
 */
void morpho_erosion(int option,
                    const VectorInt &radius,
                    const BImage& imagin,
                    BImage& imagout,
                    bool verbose)
{
  BImage imagtmp;
  int ix, iy, iz, jx, jy, jz, nx1, nx2, ny1, ny2, nz1, nz2, nbin;

  /* Initializations */

  _st_morpho_image_radius_define(radius);
  nbin = 0;
  if (verbose) nbin = morpho_count(imagin);

  /* Copy the input image into the temporary image */

  imagout = imagin;

  /* Process the structuring element along Z */

  if (RADIUS[2] != 0)
  {
    if (option == BLOCK) imagtmp = imagout;
    for (iz = 0; iz < imagin.getNDims(2); iz++)
    {
      nz1 = MIN(imagin.getNDims(2) - 1, MAX(0,iz-RADIUS[2]));
      nz2 = MIN(imagin.getNDims(2) - 1, MAX(0,iz+RADIUS[2]));
      for (iy = 0; iy < imagin.getNDims(1); iy++)
        for (ix = 0; ix < imagin.getNDims(0); ix++)
          if (! IMAGES(ix, iy, iz))
            imagout.setMaskoff(ix, iy, iz);
          else
            for (jz = nz1; jz <= nz2; jz++)
              if (!IMAGES(ix, iy, jz))
              {
                imagout.setMaskoff(ix,  iy,  iz);
                break;
              }
    }
  }

  /* Process the structuring element along Y */

  if (RADIUS[1] != 0)
  {
    if (option == BLOCK) imagtmp = imagout;
    for (iy = 0; iy < imagin.getNDims(1); iy++)
    {
      ny1 = MIN(imagin.getNDims(1) - 1, MAX(0,iy-RADIUS[1]));
      ny2 = MIN(imagin.getNDims(1) - 1, MAX(0,iy+RADIUS[1]));
      for (iz = 0; iz < imagin.getNDims(2); iz++)
        for (ix = 0; ix < imagin.getNDims(0); ix++)
          if (!IMAGES(ix, iy, iz))
            imagout.setMaskoff(ix,iy,iz);
          else
            for (jy = ny1; jy <= ny2; jy++)
              if (!IMAGES(ix, jy, iz))
              {
                imagout.setMaskoff(ix,iy,iz);
                break;
              }
    }
  }

  /* Process the structuring element along X */

  if (RADIUS[0] != 0)
  {
    if (option == BLOCK) imagtmp = imagout;
    for (ix = 0; ix < imagin.getNDims(0); ix++)
    {
      nx1 = MIN(imagin.getNDims(0) - 1, MAX(0,ix-RADIUS[0]));
      nx2 = MIN(imagin.getNDims(0) - 1, MAX(0,ix+RADIUS[0]));
      for (iz = 0; iz < imagin.getNDims(2); iz++)
        for (iy = 0; iy < imagin.getNDims(1); iy++)
          if (!IMAGES(ix, iy, iz))
            imagout.setMaskoff(ix,iy,iz);
          else
            for (jx = nx1; jx <= nx2; jx++)
              if (!IMAGES(jx, iy, iz))
              {
                imagout.setMaskoff(ix,iy,iz);
                break;
              }
    }
  }

  if (verbose) message("Erosion: %d -> %d\n", nbin, morpho_count(imagout));

  return;
}

/**
 * Perform the "dilation" of the input image
 */
void morpho_dilation(int option,
                     const VectorInt &radius,
                     const BImage& imagin,
                     BImage& imagout,
                     bool verbose)
{
  BImage imagtmp;
  int ix, iy, iz, jx, jy, jz, nx1, nx2, ny1, ny2, nz1, nz2, nbin;

  /* Initializations */

  _st_morpho_image_radius_define(radius);
  nbin = 0;
  if (verbose) nbin = morpho_count(imagin);

  /* Copy the input image into the temporary image */

  imagout = imagin;

  /* Process the structuring element along Z */

  if (RADIUS[2] != 0)
  {
    if (option == BLOCK) imagtmp = imagout;
    for (iz = 0; iz < imagin.getNDims(2); iz++)
    {
      nz1 = MIN(imagin.getNDims(2) - 1, MAX(0,iz-RADIUS[2]));
      nz2 = MIN(imagin.getNDims(2) - 1, MAX(0,iz+RADIUS[2]));
      for (iy = 0; iy < imagin.getNDims(1); iy++)
        for (ix = 0; ix < imagin.getNDims(0); ix++)
          if (IMAGES(ix, iy, iz))
            imagout.setOffset(ix,iy,iz);
          else
            for (jz = nz1; jz <= nz2; jz++)
              if (IMAGES(ix, iy, jz))
              {
                imagout.setOffset(ix,iy,iz);
                break;
              }
    }
  }

  /* Process the structuring element along Y */

  if (RADIUS[1] != 0)
  {
    if (option == BLOCK) imagtmp = imagout;
    for (iy = 0; iy < imagin.getNDims(1); iy++)
    {
      ny1 = MIN(imagin.getNDims(1) - 1, MAX(0,iy-RADIUS[1]));
      ny2 = MIN(imagin.getNDims(1) - 1, MAX(0,iy+RADIUS[1]));
      for (iz = 0; iz < imagin.getNDims(2); iz++)
        for (ix = 0; ix < imagin.getNDims(0); ix++)
          if (IMAGES(ix, iy, iz))
            imagout.setOffset(ix,iy,iz);
          else
            for (jy = ny1; jy <= ny2; jy++)
              if (IMAGES(ix, jy, iz))
              {
                imagout.setOffset(ix,iy,iz);
                break;
              }
    }
  }

  /* Process the structuring element along X */

  if (RADIUS[0] != 0)
  {
    if (option == BLOCK) imagtmp = imagout;
    for (ix = 0; ix < imagin.getNDims(0); ix++)
    {
      nx1 = MIN(imagin.getNDims(0) - 1, MAX(0,ix-RADIUS[0]));
      nx2 = MIN(imagin.getNDims(0) - 1, MAX(0,ix+RADIUS[0]));
      for (iz = 0; iz < imagin.getNDims(2); iz++)
        for (iy = 0; iy < imagin.getNDims(1); iy++)
          if (IMAGES(ix, iy, iz))
            imagout.setOffset(ix,iy,iz);
          else
            for (jx = nx1; jx <= nx2; jx++)
              if (IMAGES(jx, iy, iz))
              {
                imagout.setOffset(ix,iy,iz);
                break;
              }
    }
  }

  if (verbose) message("Dilation: %d -> %d\n", nbin, morpho_count(imagout));

  return;
}

/**
 * Perform the "intersection" of two input images
 */
void morpho_intersection(const BImage& image1,
                         const BImage& image2,
                         BImage& imagout,
                         bool verbose)
{
  int nbin1 = 0;
  int nbin2 = 0;

  if (verbose)
  {
    nbin1 = morpho_count(image1);
    nbin2 = morpho_count(image2);
  }

  for (int i = 0; i < image1.getAllocSize(); i++)
    imagout.setValue(i, image1.getValue(i) & image2.getValue(i));

  if (verbose)
    message("Intersection : %d and %d -> %d\n", nbin1, nbin2,
            morpho_count(imagout));

  return;
}

/**
 * Perform the "union" of two input images
 */
void morpho_union(const BImage& image1,
                  const BImage& image2,
                  BImage& imagout,
                  bool verbose)
{
  int nbin1 = 0;
  int nbin2 = 0;

  if (verbose)
  {
    nbin1 = morpho_count(image1);
    nbin2 = morpho_count(image2);
  }

  for (int i = 0; i < image1.getAllocSize(); i++)
    imagout.setValue(i, image1.getValue(i) | image2.getValue(i));

  if (verbose)
    message("Union: %d and %d -> %d\n", nbin1, nbin2,
            morpho_count(imagout));

  return;
}

/**
 * Perform the "negation" of the input image
 */
void morpho_negation(const BImage& imagin,
                     BImage& imagout,
                     bool verbose)
{
  int nbin = 0;

  if (verbose) nbin = morpho_count(imagin);

  for (int i = 0; i < imagin.getAllocSize(); i++)
    imagout.setValue(i, ~imagin.getValue(i));

  if (verbose) message("Negation: %d -> %d\n", nbin, morpho_count(imagout));

  return;
}

/**
 * Returns the volume (number of pixels) assigned to the grain (Binary Image == 1)
 */
int morpho_count(const BImage& imagin)
{
  int ncount = 0;
  for (int iz = 0; iz < imagin.getNDims(2); iz++)
    for (int iy = 0; iy < imagin.getNDims(1); iy++)
      for (int ix = 0; ix < imagin.getNDims(0); ix++)
        if (imagin.getValue(ix, iy, iz)) ncount++;
  return (ncount);
}

/**
 * Perform the "opening" of the input image (i.e. an erosion followed by a dilation)
 */
void morpho_opening(int option,
                    const VectorInt &radius,
                    const BImage& imagin,
                    BImage& imagout,
                    bool verbose)
{
  BImage imagtmp = imagin;

  morpho_erosion(option, radius, imagin, imagtmp, verbose);

  morpho_dilation(option, radius, imagtmp, imagout, verbose);

  return;
}

/**
 * Perform the "closing" of the input image (i.e. a dilation followed by an erosion)
 */
void morpho_closing(int option,
                    const VectorInt &radius,
                    const BImage& imagin,
                    BImage& imagout,
                    bool verbose)
{
  BImage imagtmp = imagin;

  morpho_dilation(option, radius, imagin, imagtmp, verbose);

  morpho_erosion(option, radius, imagtmp, imagout, verbose);

  return;
}

/**
 * Converts an input image (double) into an image (in place)
 */
void morpho_double2imageInPlace(const VectorInt &nx,
                                const VectorDouble &tab,
                                double vmin,
                                double vmax,
                                BImage &imagout,
                                bool verbose)
{
  imagout.init(nx);
  VectorInt NX = imagout.getNDimsExt(3);
  unsigned char mot = 0;
  int ind = 0;
  int ecr = 0;
  for (int iz = 0; iz < imagout.getNDims(2); iz++)
    for (int iy = 0; iy < imagout.getNDims(1); iy++)
      for (int ix = 0; ix < imagout.getNDims(0); ix++)
      {
        ind++;
        double val = TAB(ix, iy, iz);
        double result = 1.;
        if (FFFF(val)) result = 0.;
        if (!FFFF(vmin) && val < vmin) result = 0.;
        if (!FFFF(vmax) && val >= vmax) result = 0.;
        mot = (mot << 1) + (unsigned char) result;
        if (ind == 8)
        {
          imagout.setValue(ecr, mot);
          mot = 0;
          ind = 0;
          ecr++;
        }
      }

  if (ind != 0) imagout.setValue(ecr, (mot << (8 - ind)));

  if (verbose)
    message("Translation: %d  / %d\n", morpho_count(imagout), imagout.getNPixels());

  return;
}

/**
 * Converts an input image (double) into a returned image
 */
BImage morpho_double2image(const VectorInt &nx,
                           const VectorDouble &tab,
                           double vmin,
                           double vmax,
                           bool verbose)
{
  BImage imagout(nx);
  morpho_double2imageInPlace(nx, tab, vmin, vmax, imagout, verbose);
  return imagout;
}

/**
 * Converts an image into an array (double)
 */
void morpho_image2double(const BImage& imagin,
                         int mode,
                         double grain,
                         double pore,
                         VectorDouble &tab,
                         bool verbose)
{
  VectorInt NX = imagin.getNDimsExt(3);
  int nxyz = imagin.getNPixels();
  if (verbose)
    message("Translation: %d / %d\n", morpho_count(imagin), nxyz);

  for (int iz = 0; iz < imagin.getNDims(2); iz++)
    for (int iy = 0; iy < imagin.getNDims(1); iy++)
      for (int ix = 0; ix < imagin.getNDims(0); ix++)
      {
        double value = imagin.getValue(ix,iy,iz) ? grain : pore;
        switch (mode)
        {
          case 0:
            TAB(ix,iy,iz) = value;
            break;

          case 1:
            TAB(ix,iy,iz) += value;
            break;

          case -1:
            TAB(ix,iy,iz) -= value;
            break;
        }
      }
  return;
}

/**
 * Returns the vector of distances from the grain to the edge
 */
void morpho_distance(int option,
                     const VectorInt &radius,
                     bool flagDistErode,
                     BImage& imagin,
                     VectorDouble &dist,
                     bool verbose)
{
  BImage imagout = imagin;
  int nxyz = imagin.getNPixels();

  /* Copy the initial image in the distance array */

  morpho_image2double(imagin, 0, 1, 0, dist);

  /* Processing loop */

  int iter = 0;
  if (flagDistErode)
  {
    while (morpho_count(imagin) != 0)
    {
      iter++;
      morpho_erosion(option, radius, imagin, imagout);
      morpho_duplicate(imagout, imagin);
      morpho_image2double(imagin, 1, 1, 0, dist);
      if (verbose)
      {
        int count = morpho_count(imagin);
        message("Iteration %d: Current (%d/%d)\n",iter,count,nxyz);
      }
    }
  }
  else
  {
    int incr = 0;
    while (morpho_count(imagin) != nxyz)
    {
      iter++;
      incr++;
      morpho_dilation(option, radius, imagin, imagout);
      morpho_duplicate(imagout, imagin);
      morpho_image2double(imagin, -1, 1, 0, dist);
      if (verbose)
      {
        int count = morpho_count(imagin);
        message("Iteration %d: Current (%d/%d)\n",iter,count,nxyz);
      }
    }
    morpho_image2double(imagin, 1, incr, 0, dist);
  }

  /* Turn the distance to a positive value */

  for (int i = 0; i < nxyz; i++)
    dist[i] = ABS(dist[i]);

  return;
}

/**
 * Create the array of index shifts (relatively to the center cell) for a dilation by 'radius' of a
 * regular grid
 */
VectorInt gridcell_neigh(int ndim,
                         int option,
                         int radius,
                         bool flag_center,
                         bool verbose)
{
  int *indg0, *indg1, ecr, flag_count, nech;
  VectorInt indret;

  /* Initializations */

  int nvois = 0;
  indg0 = indg1 = nullptr;

  /* Create the grid attributes */

  VectorInt    nx(ndim);
  VectorDouble x0(ndim);
  VectorDouble dx(ndim);
  for (int idim = 0; idim < ndim; idim++)
  {
    dx[idim] = 1.;
    nx[idim] = 1 + 2 * radius;
    x0[idim] = 0.;
  }

  DbGrid* grid = DbGrid::create(nx, dx, x0, VectorDouble(), ELoadBy::SAMPLE,
                                VectorDouble(), VectorString(),
                                VectorString(), 1);

  /* (Maximum) core allocation */

  nech = grid->getSampleNumber();
  indret.resize(nech * ndim);
  indg0 = db_indg_alloc(grid);
  if (indg0 == nullptr) goto label_end;
  indg1 = db_indg_alloc(grid);
  if (indg1 == nullptr) goto label_end;

  /* Scan the grid nodes */

  ecr = 0;
  db_index_sample_to_grid(grid, nech / 2, indg0);
  for (int iech = 0; iech < nech; iech++)
  {
    db_index_sample_to_grid(grid, iech, indg1);
    flag_count = 0;
    for (int idim = 0; idim < ndim; idim++)
    {
      indg1[idim] -= indg0[idim];
      if (indg1[idim] != 0) flag_count++;
    }

    if (flag_center && iech == nech / 2) continue;
    if (option == CROSS && flag_count > 1) continue;

    /* Add this cell */

    for (int idim = 0; idim < ndim; idim++)
      indret[ecr++] = indg1[idim];
  }

  /* Resizing the returned array */

  indret.resize(ecr);
  nvois = ecr / ndim;

  /* Optional printout */

  if (verbose && nvois > 0)
  {
    ecr = 0;
    message("Grid Dilation: %d samples\n", nvois);
    for (int i = 0; i < nvois; i++)
    {
      message("  Neigh %3d:", i + 1);
      for (int idim = 0; idim < ndim; idim++)
        message(" %2d", indret[ecr++]);
      message("\n");
    }
  }

  label_end:
  if (grid != nullptr) delete grid;
  indg0 = db_indg_free(indg0);
  indg1 = db_indg_free(indg1);
  return (indret);
}

/**
 * Calculate the gradient orientations of a colored image
 */
void _morpho_angle2D(DbGrid *dbgrid, const VectorInt &radius, int iptr0)
{
  int iad;
  int pivot = 0;
  double a[3], b[2], x[2], result;

  /* Initializations */

  VectorInt NX = dbgrid->getNXsExt(2);
  _st_morpho_image_radius_define(radius);
  int ndim = dbgrid->getNDim();
  VectorInt indg(ndim);
  int iptrz = dbgrid->getColIdxByLocator(ELoc::Z);

  /* Processing */

  int iiz = 0;
  for (int iiy = 0; iiy < NX[1]; iiy++)
    for (int iix = 0; iix < NX[0]; iix++)
    {
      for (int i = 0; i < 3; i++)  a[i] = 0.;
      for (int i = 0; i < 2; i++)
        b[i] = x[i] = 0.;
      indg[0] = iix;
      indg[1] = iiy;
      iad = dbgrid->indiceToRank(indg);
      double z0 = dbgrid->getArray(iad, iptrz);
      if (FFFF(z0)) continue;

      for (int jx = -RADIUS[0]; jx <= RADIUS[0]; jx++)
        for (int jy = -RADIUS[1]; jy <= RADIUS[1]; jy++)
        {
          int ix = iix + jx;
          if (ix < 0 || ix >= NX[0]) continue;
          int iy = iiy + jy;
          if (iy < 0 || iy >= NX[1]) continue;
          indg[0] = ix;
          indg[1] = iy;
          iad = dbgrid->indiceToRank(indg);
          double xi = dbgrid->getCoordinate(iad, 0);
          double yi = dbgrid->getCoordinate(iad, 1);
          double zi = dbgrid->getArray(iad, iptrz);
          if (FFFF(zi)) continue;
          zi -= z0;
          a[0] += xi * xi;
          a[1] += xi * yi;
          a[2] += yi * yi;
          b[0] += xi * zi;
          b[1] += yi * zi;
        }
      if (matrix_solve(0, a, b, x, 2, 1, &pivot)) continue;

      result = ut_rad2deg(atan2(x[1], x[0]));
      result += 90.;
      while (result < -180)
        result += 360.;
      while (result > 180)
        result -= 360.;

      if (ndim >= 1) indg[0] = iix;
      if (ndim >= 2) indg[1] = iiy;
      if (ndim >= 3) indg[2] = iiz;
      dbgrid->setArray(dbgrid->indiceToRank(indg), iptr0, result);
    }
  return;
}

/**
 * Calculate the gradient components of a colord image
 */
void _morpho_gradients(DbGrid *dbgrid, int iptr)
{
  int j1, j2, number;

  /* Preliminary check */

  VectorInt NX = dbgrid->getNXsExt(3);
  int ndim = dbgrid->getNDim();
  VectorInt indg(ndim);
  int iptrz = dbgrid->getColIdxByLocator(ELoc::Z);

  /* Calculate the Gradient components */

  for (int ix = 0; ix < NX[0]; ix++)
    for (int iy = 0; iy < NX[1]; iy++)
      for (int iz = 0; iz < NX[2]; iz++)
      {
        for (int idim = 0; idim < ndim; idim++)
        {
          if (ndim >= 1) indg[0] = ix;
          if (ndim >= 2) indg[1] = iy;
          if (ndim >= 3) indg[2] = iz;

          int nmax = dbgrid->getNX(idim);
          double dinc = dbgrid->getDX(idim);

          double v1 = 0.;
          double v2 = 0.;
          number = 0;
          if (idim == 0)
          {
            j1 = (ix + 1 > nmax - 1) ? ix : ix + 1;
            indg[0] = j1;
            v1 = dbgrid->getArray(dbgrid->indiceToRank(indg), iptrz);
            if (FFFF(v1)) continue;
            j2 = (ix - 1 < 0) ? ix : ix - 1;
            indg[0] = j2;
            v2 = dbgrid->getArray(dbgrid->indiceToRank(indg), iptrz);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          if (idim == 1)
          {
            j1 = (iy + 1 > nmax - 1) ? iy : iy + 1;
            indg[1] = j1;
            v1 = dbgrid->getArray(dbgrid->indiceToRank(indg), iptrz);
            if (FFFF(v1)) continue;
            j2 = (iy - 1 < 0) ? iy : iy - 1;
            indg[1] = j2;
            v2 = dbgrid->getArray(dbgrid->indiceToRank(indg), iptrz);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          if (idim == 2)
          {
            j1 = (iz + 1 > nmax - 1) ? iz : iz + 1;
            indg[2] = j1;
            v1 = dbgrid->getArray(dbgrid->indiceToRank(indg), iptrz);
            if (FFFF(v1)) continue;
            j2 = (iz - 1 < 0) ? iz : iz - 1;
            indg[2] = j2;
            v2 = dbgrid->getArray(dbgrid->indiceToRank(indg), iptrz);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          if (ndim >= 1) indg[0] = ix;
          if (ndim >= 2) indg[1] = iy;
          if (ndim >= 3) indg[2] = iz;
          double delta = (v1 - v2) / (number * dinc);
          dbgrid->setArray(dbgrid->indiceToRank(indg), iptr + idim, delta);
        }
      }
}

/**
 * Perform a morphological operation with a DbGrid
 */
int _db_morpho_calc(DbGrid *dbgrid,
                    int iptr0,
                    const EMorpho &oper,
                    double vmin,
                    double vmax,
                    int option,
                    const VectorInt &radius,
                    bool flagDistErode,
                    bool verbose)
{
  int ntotal = dbgrid->getSampleNumber();
  VectorInt nxy = dbgrid->getNXs();

  VectorDouble tabin = dbgrid->getColumnByLocator(ELoc::Z);
  BImage image = morpho_double2image(nxy,tabin,vmin,vmax);

  if (verbose)
  {
    message("Morphological operation = %s\n",oper.getDescr().c_str());
    message("Initial image = %d/%d\n",morpho_count(image),ntotal);
  }

  bool alreadyLoaded = false;
  bool alreadySaved = false;
  BImage image2 = BImage(nxy);
  VectorDouble tabout = VectorDouble(ntotal, TEST);
  if (oper == EMorpho::THRESH)
  {
    morpho_duplicate(image, image2);
  }
  else if (oper == EMorpho::NEGATION)
  {
    morpho_negation(image, image2);
  }
  else if (oper == EMorpho::EROSION)
  {
    morpho_erosion(option, radius, image, image2);
  }
  else if (oper == EMorpho::DILATION)
  {
    morpho_dilation(option, radius, image, image2);
  }
  else if (oper == EMorpho::OPEN)
  {
    morpho_opening(option, radius, image, image2);
  }
  else if (oper == EMorpho::CLOSE)
  {
    morpho_closing(option, radius, image, image2);
  }
  else if (oper == EMorpho::CC)
  {
    tabout = morpho_labelling(0, 0, image, TEST, verbose);
    alreadyLoaded = true;
  }
  else if (oper == EMorpho::CCSIZE)
  {
    tabout = morpho_labelling(0, 1, image, TEST, verbose);
    alreadyLoaded = true;
  }
  else if (oper == EMorpho::DISTANCE)
  {
    morpho_distance(option,radius,flagDistErode,image,tabout,verbose);
    alreadyLoaded = true;
  }
  else if (oper == EMorpho::ANGLE)
  {
    _morpho_angle2D(dbgrid,radius,iptr0);
    alreadyLoaded = true;
    alreadySaved = true;
  }
  else if (oper == EMorpho::GRADIENT)
  {
    _morpho_gradients(dbgrid, iptr0);
    alreadyLoaded = true;
    alreadySaved = true;
  }
  else
  {
    messerr("Not programmed yet\n");
    return 1;
  }

  if (! alreadyLoaded)
  {
    if (verbose)
      message("Resulting image = %d/%d\n",morpho_count(image2),ntotal);
    morpho_image2double(image2, 0, 1., 0., tabout);
  }
  if (! alreadySaved)
  {
    dbgrid->setColumnByUID(tabout, iptr0);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Evaluates the spill point
 **
 ** \return  The Spill_Res structure which contains:
 ** \return   h      elevation of the spill point
 ** \return   th     maximum reservoir thickness
 ** \return   ix0    location of the spill point grid node along X
 ** \return   iy0    location of the spill point grid node along Y
 **
 ** \param[in]  dbgrid        Grid Db structure
 ** \param[in]  name_depth    Name of the variable containing the depth
 ** \param[in]  name_data     Name of the variable containing the data
 ** \param[in]  option        0 for 4-connectivity; 1 for 8-connectivity
 ** \param[in]  flag_up       TRUE when working in elevation; FALSE in depth
 ** \param[in]  verbose_step  Step for the verbose flag
 ** \param[in]  hmax          maximum reservoir thickness (FFFF not used)
 **
 ** \remark  The variable 'ind_data', which contains the constraints, must
 ** \remark  be set to:
 ** \remark  0 for an idle node
 ** \remark  1 for a node located outside the reservoir
 ** \remark  2 for a node belonging to the reservoir
 ** \remark  The numbering of the grid node corresponding to the spill point
 ** \remark  must start with 1
 **
 *****************************************************************************/
Spill_Res spillPoint(DbGrid *dbgrid,
                     const String& name_depth,
                     const String& name_data,
                     int option,
                     bool flag_up,
                     int verbose_step,
                     double hmax)
{
  Spill_Res res;
  double h, th;
  int ix0, iy0;

  int ind_depth = dbgrid->getUID(name_depth);
  int ind_data = dbgrid->getUID(name_data);
  if (ind_depth < 0 || ind_data < 0)
  {
    messerr("Variables 'name_depth' and 'name_data' are compulsory");
    res.success = false;
    return res;
  }

  int error = spill_point(dbgrid, ind_depth, ind_data, option, flag_up,
                          verbose_step, hmax, &h, &th, &ix0, &iy0);

  res.success = (error) ? false : true;
  res.h = h;
  res.th = th;
  res.ix0 = ix0;
  res.iy0 = iy0;

  return res;
}
