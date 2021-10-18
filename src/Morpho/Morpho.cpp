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
#include "Morpho/Morpho.hpp"
#include "Basic/Utilities.hpp"
#include "Db/Db.hpp"
#include "geoslib_f.h"
#include "geoslib_old_f.h"

static int NX[3],NXYZ,NRED,RADIUS[3];
static int LARGE = 9999999;
static unsigned char Offset[]  = { 128,  64,  32,  16,   8,   4,   2,   1 };
static unsigned char Maskoff[] = { 127, 191, 223, 239, 247, 251, 253, 254 };

/*! \cond */
#define CROSS 0
#define BLOCK 1
#define GRID_ADD(ix,iy,iz)  ((ix)+(NX[0]*((iy)+NX[1]*(iz))))
#define IN_GRID(ix,iy,iz)    (ix>=0 && ix<NX[0] &&  \
                              iy>=0 && iy<NX[1] &&  \
                              iz>=0 && iz<NX[2])
#define DIVIDE(i,j,k)       (GRID_ADD(i,j,k) / 8)
#define RESIDU(i,j,k)       (GRID_ADD(i,j,k) % 8)
#define OFFSET(i,j,k)       (Offset [RESIDU(i,j,k)])
#define MASKOFF(i,j,k)      (Maskoff[RESIDU(i,j,k)])
#define IMAGIN(i,j,k)       (imagin [DIVIDE(i,j,k)])
#define IMAGOUT(i,j,k)      (imagout[DIVIDE(i,j,k)])
#define IMAGTMP(i,j,k)      (imagtmp[DIVIDE(i,j,k)])
#define VALIN(ix,iy,iz)     (IMAGIN (ix,iy,iz) & OFFSET(ix,iy,iz))
#define VALTMP(ix,iy,iz)    (IMAGTMP(ix,iy,iz) & OFFSET(ix,iy,iz))
#define IMAGES(ix,iy,iz)    ((option == BLOCK) ? VALTMP(ix,iy,iz) : VALIN(ix,iy,iz))

#define TAB(i,j,k)          (tab[GRID_ADD(i,j,k)])
#define TABOUT(i,j,k)       (tabout[GRID_ADD(i,j,k)])
#define COMPNUM(i,j,k)      (compnum[GRID_ADD(i,j,k)])
/*! \endcond */

/*****************************************************************************/
/*!
**  Defines the image size (global variables)
**
** \param[in]  nx     Number of grid meshes
**
*****************************************************************************/
void _st_morpho_image_size_define(const VectorInt& nx)
{
  int size;

  size = static_cast<int> (nx.size());
  NX[0] = (size > 0) ? nx[0] : 1;
  NX[1] = (size > 1) ? nx[1] : 1;
  NX[2] = (size > 2) ? nx[2] : 1;

  NXYZ  = NX[0] * NX[1] * NX[2];
  NRED  = ((NXYZ - 8) / 8 + 1);
}

/*****************************************************************************/
/*!
**  Defines the image radius (global variables)
**
** \param[in]  radius Radius of the structural element
**
*****************************************************************************/
void _st_morpho_image_radius_define(const VectorInt& radius)
{
  int size = static_cast<int> (radius.size());
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
int _st_morpho_label_size(const VectorDouble& compnum,
                          int nbcomp,
                          VectorInt& sizes)
{
  int total = 0;
  for (int i=0; i<NXYZ; i++)
  {
    int val = (int) compnum[i];
    if (val > 0 && val <= nbcomp)
    {
      sizes[val-1]++;
      total++;
    }
  }
  return(total);
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
void _st_morpho_label_order(VectorDouble& compnum,
                            const VectorInt& order,
                            int nbcomp)
{
  int i,j,val,found;

  for (i=0; i<NXYZ; i++)
  {
    val = (int) compnum[i];
    if (val <= 0) continue;
    for (j=nbcomp-1,found= -1; j>=0 && found<0; j--)
      if (val == order[j]) found = j;
    if (found < 0) messageAbort("st_morpho_label_order");
    compnum[i] = nbcomp - found;
  }
  return;
}

/*****************************************************************************/
/*!
**  Returns the dimension of an image
**
** \return  Dimension of an image (integer)
**
** \param[out]  nx     Number of grid meshes (dimension = 3)
**
*****************************************************************************/
int morpho_image_size(const VectorInt& nx)

{
  _st_morpho_image_size_define(nx);
  return(NRED);
}

/*****************************************************************************/
/*!
**  Copy the input image into the output image
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  imagin  input image
**
** \param[out] imagout output image
**
*****************************************************************************/
void morpho_duplicate(const VectorInt& nx,
                      const VectorUChar& imagin,
                      VectorUChar& imagout)
{
  int i;

  _st_morpho_image_size_define(nx);
  for (i=0; i<NRED; i++) imagout[i] = imagin[i];

  return;
}

/*****************************************************************************/
/*!
**  Labels the connected components for a 3D image
**
** \return  Number of connected components
**
** \param[in]  nx        Number of grid meshes (dimension = 3)
** \param[in]  option    connectivity option (CROSS or BLOCK)
** \param[in]  flag_size 1 if cell contains volume of connected component
**                       0 if cell contains rank of connected component
** \param[in]  imagin    input image
** \param[in]  ccvoid    Value assigned to pixels which do not belong to any
**                       connected component
** \param[in]  verbose   Verbose flag
**
** \param[out] compnum   output array containing the component index
**
** \remark  The labels are sorted by decreasing sizes and an optional message
** \remark  is issued for displaying the component sizes
**
*****************************************************************************/
int morpho_labelling(const VectorInt& nx,
                     int option,
                     int flag_size,
                     const VectorUChar& imagin,
                     double ccvoid,
                     VectorDouble& compnum,
                     bool verbose)
{
  int  jx,jy,jz,total,count,ref,iad,local,part_grain;
  int  i,size,nbtest,nbcomp,ix,iy,iz,icomp,ival,itest,il,jcomp[26];
  VectorInt list_array,sizes,order;
  int id[26][3] = {
    {-1,00,00}, {00,-1,00}, {00,00,-1}, {01,00,00}, {00,01,00}, {00,00,01},
    {-1,-1,-1}, {00,-1,-1}, {01,-1,-1}, {-1,00,-1}, {-1,01,-1}, {00,01,-1},
    {01,00,-1}, {01,01,-1},
    {-1,-1,00}, {-1,01,00}, {01,-1,00}, {01,01,00},
    {-1,-1,01}, {00,-1,01}, {01,-1,01}, {-1,00,01}, {-1,01,01}, {00,01,01},
    {01,00,01}, {01,01,01} };
  int ndel[2]  = { 6, 26 };
  int quantum = 100;

  /* Initializations */

  _st_morpho_image_size_define(nx);
  nbtest = ndel[option];
  nbcomp = total = 0;

  /* Attempt to allocate the initial quantum */

  size = MIN(quantum,NXYZ);
  list_array.resize(size);
  for (i=0; i<NXYZ; i++) compnum[i] = 0.;

  /* Numbering the connected components */

  for (iz=0; iz<NX[2]; iz++)
    for (iy=0; iy<NX[1]; iy++)
      for (ix=0; ix<NX[0]; ix++)
      {
        if (VALIN(ix,iy,iz))
        {
          icomp = LARGE;
          for (itest=0; itest<nbtest; itest++)
          {
            jcomp[itest] = LARGE;
            jx = ix + id[itest][0];
            jy = iy + id[itest][1];
            jz = iz + id[itest][2];
            if (IN_GRID(jx,jy,jz))
            {
              ival = (int) COMPNUM(jx,jy,jz);
              if (ival > 0)
              {
                jcomp[itest] = list_array[ival-1];
                if (jcomp[itest] < icomp) icomp = jcomp[itest];
              }
            }
          }
          if (icomp != LARGE)
          {
            for (itest=0; itest<nbtest; itest++)
              if (jcomp[itest] != LARGE && jcomp[itest] != icomp)
                for (il=jcomp[itest]; il<=nbcomp; il++)
                  if (list_array[il-1] == jcomp[itest])
                    list_array[il-1] = icomp;
          }
          else
          {
            nbcomp++;
            if (nbcomp > size)
            {
              size += quantum;
              list_array.resize(size);
            }
            list_array[nbcomp-1] = nbcomp;
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

  for (il=ival=0; il<nbcomp; il++)
    if( list_array[il] == il+1)
      list_array[il] = ++ival;
    else
      list_array[il] = list_array[list_array[il]-1];
  nbcomp = ival;

  /* Update the planes */

  for (iz=0; iz<NX[2]; iz++)
    for (iy=0; iy<NX[1]; iy++)
      for (ix=0; ix<NX[0]; ix++)
      {
        iad = (int) COMPNUM(ix,iy,iz);
        if (iad != 0)
          COMPNUM(ix,iy,iz) = list_array[iad-1];
        else
          COMPNUM(ix,iy,iz) = TEST;
      }

  /* Order the components */

  if (nbcomp > 0)
  {
    sizes.resize(nbcomp,0);
    total = _st_morpho_label_size(compnum,nbcomp,sizes);
    order.resize(nbcomp);
    for (i=0; i<nbcomp; i++) order[i] = i+1;
    ut_sort_int(1,nbcomp,order.data(),sizes.data());
    _st_morpho_label_order(compnum,order,nbcomp);

    if (flag_size)
      for (i=0; i<NXYZ; i++)
      {
        if (FFFF(compnum[i])) continue;
        ival = (int) compnum[i];
        compnum[i] = sizes[ival-1];
      }

    for (i=0; i<NXYZ; i++)
      if (FFFF(compnum[i])) compnum[i] = ccvoid;
  }

  /* Display the label component size statistics */

  if (verbose)
  {
    message("Labelling: %d\n",total);
    if (nbcomp == 0)
    {
      message("No grain has been detected\n");
    }
    else
    {
      message("Number of connected components = %d\n\n", nbcomp);
      message("   Component     Number         Total     Cumul (percent)\n");
      ref = sizes[order[nbcomp-1]-1];
      count = part_grain = 0;
      for (i=nbcomp-1; i>=0; i--)
      {
        local = sizes[order[i]-1];

        if (local != ref)
        {
          message("%12d  %9d  %12d      %7.3f\n",
                  ref, count, count*ref, 100.*part_grain/total);
          ref = local;
          count = 0;
        }

        part_grain += local;
        count++;
      }
      message("%12d  %9d  %12d      %7.3f\n",
              ref, count, count*ref, 100.*part_grain/total);
    }
  }

  return nbcomp;
}

/*****************************************************************************/
/*!
**  Labels the connected components for a 3D image
**
** \return  The array of connected components newly created
**
** \param[in]  nx        Number of grid meshes (dimension = 3)
** \param[in]  option    connectivity option (CROSS or BLOCK)
** \param[in]  flag_size 1 if cell contains volume of connected component
**                       0 if cell contains rank of connected component
** \param[in]  imagin    input image
** \param[in]  ccvoid    Value assigned to pixels which do not belong to any
**                       connected component
** \param[in]  verbose   Verbose flag
**
** \remark  The labels are sorted by decreasing sizes and an optional message
** \remark  is issued for displaying the component sizes
**
*****************************************************************************/
VectorDouble morpho_labelling(const VectorInt& nx,
                              int option,
                              int flag_size,
                              const VectorUChar& imagin,
                              double ccvoid,
                              bool verbose)
{
  _st_morpho_image_size_define(nx);
  VectorDouble compnum(NXYZ);
  (void) morpho_labelling(nx,option,flag_size,imagin,ccvoid,compnum,verbose);
  return compnum;
}

/*****************************************************************************/
/*!
**  Returns the array of dimensions of the connex components
**
** \return  Vector of sizes
**
** \param[in]  nx        Number of grid meshes (dimension = 3)
** \param[in]  option    connectivity option (CROSS or BLOCK)
** \param[in]  imagin    input image
**
** \remark  The labels are sorted by decreasing sizes and an optional message
** \remark  is issued for displaying the component sizes
**
*****************************************************************************/
VectorInt morpho_labelsize(const VectorInt& nx, int option, const VectorUChar& imagin)
{
  VectorInt sizes;

  _st_morpho_image_size_define(nx);
  VectorDouble compnum(NXYZ);
  int nbcomp = morpho_labelling(nx,option,0,imagin,TEST,compnum);
  if (nbcomp > 0)
  {
    sizes.resize(nbcomp,0);
    (void) _st_morpho_label_size(compnum,nbcomp,sizes);
  }
  return sizes;
}

/*****************************************************************************/
/*!
**  Performs a morphological erosion
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  option  Option of the structuring element (CROSS pr BLOCK)
** \param[in]  radius  Radius of the structuring element (dimension = 3)
** \param[in]  imagin  input image
** \param[in]  verbose Verbose flag
**
** \param[out] imagout Output image
**
*****************************************************************************/
void morpho_erosion(const VectorInt& nx,
                    int option,
                    const VectorInt& radius,
                    const VectorUChar& imagin,
                    VectorUChar& imagout,
                    bool verbose)
{
  VectorUChar imagtmp;
  int    ix,iy,iz,jx,jy,jz,nx1,nx2,ny1,ny2,nz1,nz2,nbin;

  /* Initializations */

  _st_morpho_image_size_define(nx);
  _st_morpho_image_radius_define(radius);
  nbin = 0;
  if (verbose) nbin = morpho_count(nx,imagin);

  /* Copy the input image into the temporary image */

  if (option == BLOCK) imagtmp = morpho_image_manage(nx);
  morpho_duplicate(nx,imagin,imagout);

  /* Process the structuring element along Z */

  if (RADIUS[2] != 0)
  {
    if (option == BLOCK) morpho_duplicate(nx,imagout,imagtmp);
    for (iz=0; iz<NX[2]; iz++)
    {
      nz1 = MIN(NX[2]-1,MAX(0,iz-RADIUS[2]));
      nz2 = MIN(NX[2]-1,MAX(0,iz+RADIUS[2]));
      for (iy=0; iy<NX[1]; iy++)
        for (ix=0; ix<NX[0]; ix++)
          if (! IMAGES(ix,iy,iz))
            IMAGOUT(ix,iy,iz) &= MASKOFF(ix,iy,iz);
          else
            for (jz=nz1; jz<=nz2; jz++)
              if (! IMAGES(ix,iy,jz))
              {
                IMAGOUT(ix,iy,iz) &= MASKOFF(ix,iy,iz);
                break;
              }
    }
  }

  /* Process the structuring element along Y */

  if (RADIUS[1] != 0)
  {
    if (option == BLOCK) morpho_duplicate(nx,imagout,imagtmp);
    for (iy=0; iy<NX[1]; iy++)
    {
      ny1 = MIN(NX[1]-1,MAX(0,iy-RADIUS[1]));
      ny2 = MIN(NX[1]-1,MAX(0,iy+RADIUS[1]));
      for (iz=0; iz<NX[2]; iz++)
        for (ix=0; ix<NX[0]; ix++)
          if (! IMAGES(ix,iy,iz))
            IMAGOUT(ix,iy,iz) &= MASKOFF(ix,iy,iz);
          else
            for (jy=ny1; jy<=ny2; jy++)
              if (! IMAGES(ix,jy,iz))
              {
                IMAGOUT(ix,iy,iz) &= MASKOFF(ix,iy,iz);
                break;
              }
    }
  }

  /* Process the structuring element along X */

  if (RADIUS[0] != 0)
  {
    if (option == BLOCK) morpho_duplicate(nx,imagout,imagtmp);
    for (ix=0; ix<NX[0]; ix++)
    {
      nx1 = MIN(NX[0]-1,MAX(0,ix-RADIUS[0]));
      nx2 = MIN(NX[0]-1,MAX(0,ix+RADIUS[0]));
      for (iz=0; iz<NX[2]; iz++)
        for (iy=0; iy<NX[1]; iy++)
          if (! IMAGES(ix,iy,iz))
            IMAGOUT(ix,iy,iz) &= MASKOFF(ix,iy,iz);
          else
            for (jx=nx1; jx<=nx2; jx++)
              if (! IMAGES(jx,iy,iz))
              {
                IMAGOUT(ix,iy,iz) &= MASKOFF(ix,iy,iz);
                break;
              }
    }
  }

  if (verbose)
    message("Erosion: %d -> %d\n",nbin,morpho_count(nx,imagout));

  return;
}

/*****************************************************************************/
/*!
**  Performs a morphological dilation
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  option  Option of the structuring element (CROSS pr BLOCK)
** \param[in]  radius  Radius of the structuring element (dimension = 3)
** \param[in]  imagin  input image
** \param[in]  verbose Verbose flag
**
** \param[out] imagout output image
**
*****************************************************************************/
void morpho_dilation(const VectorInt& nx,
                     int option,
                     const VectorInt& radius,
                     const VectorUChar& imagin,
                     VectorUChar& imagout,
                     bool verbose)
{
  VectorUChar imagtmp;
  int    ix,iy,iz,jx,jy,jz,nx1,nx2,ny1,ny2,nz1,nz2,nbin;

  /* Initializations */

  _st_morpho_image_size_define(nx);
  _st_morpho_image_radius_define(radius);
  nbin = 0;
  if (verbose) nbin = morpho_count(nx,imagin);

  /* Copy the input image into the temporary image */

  if (option == BLOCK) imagtmp = morpho_image_manage(nx);
  morpho_duplicate(nx,imagin,imagout);

  /* Process the structuring element along Z */

  if (RADIUS[2] != 0)
  {
    if (option == BLOCK) morpho_duplicate(nx,imagout,imagtmp);
    for (iz=0; iz<NX[2]; iz++)
    {
      nz1 = MIN(NX[2]-1,MAX(0,iz-RADIUS[2]));
      nz2 = MIN(NX[2]-1,MAX(0,iz+RADIUS[2]));
      for (iy=0; iy<NX[1]; iy++)
        for (ix=0; ix<NX[0]; ix++)
          if (IMAGES(ix,iy,iz))
            IMAGOUT(ix,iy,iz) |= OFFSET(ix,iy,iz);
          else
            for (jz=nz1; jz<=nz2; jz++)
              if (IMAGES(ix,iy,jz))
              {
                IMAGOUT(ix,iy,iz) |= OFFSET(ix,iy,iz);
                break;
              }
    }
  }

  /* Process the structuring element along Y */

  if (RADIUS[1] != 0)
  {
    if (option == BLOCK) morpho_duplicate(nx,imagout,imagtmp);
    for (iy=0; iy<NX[1]; iy++)
    {
      ny1 = MIN(NX[1]-1,MAX(0,iy-RADIUS[1]));
      ny2 = MIN(NX[1]-1,MAX(0,iy+RADIUS[1]));
      for (iz=0; iz<NX[2]; iz++)
        for (ix=0; ix<NX[0]; ix++)
          if (IMAGES(ix,iy,iz))
            IMAGOUT(ix,iy,iz) |= OFFSET(ix,iy,iz);
          else
            for (jy=ny1; jy<=ny2; jy++)
              if (IMAGES(ix,jy,iz))
              {
                IMAGOUT(ix,iy,iz) |= OFFSET(ix,iy,iz);
                break;
              }
    }
  }

  /* Process the structuring element along X */

  if (RADIUS[0] != 0)
  {
    if (option == BLOCK) morpho_duplicate(nx,imagout,imagtmp);
    for (ix=0; ix<NX[0]; ix++)
    {
      nx1 = MIN(NX[0]-1,MAX(0,ix-RADIUS[0]));
      nx2 = MIN(NX[0]-1,MAX(0,ix+RADIUS[0]));
      for (iz=0; iz<NX[2]; iz++)
        for (iy=0; iy<NX[1]; iy++)
          if (IMAGES(ix,iy,iz))
            IMAGOUT(ix,iy,iz) |= OFFSET(ix,iy,iz);
          else
            for (jx=nx1; jx<=nx2; jx++)
              if (IMAGES(jx,iy,iz))
              {
                IMAGOUT(ix,iy,iz) |= OFFSET(ix,iy,iz);
                break;
              }
    }
  }

  if (verbose)
    message("Dilation: %d -> %d\n",nbin,morpho_count(nx,imagout));

  return;
}

/*****************************************************************************/
/*!
**  Performs the intersection of two binary images
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  image1  first input image
** \param[in]  image2  second input image
** \param[in]  verbose Verbose flag
**
** \param[out] imagout output image
**
*****************************************************************************/
void morpho_intersection(const VectorInt& nx,
                         const VectorUChar& image1,
                         const VectorUChar& image2,
                         VectorUChar& imagout,
                         bool verbose)
{
  int i,nbin1,nbin2;

  nbin1 = nbin2 = 0;
  _st_morpho_image_size_define(nx);

  if (verbose)
  {
    nbin1 = morpho_count(nx,image1);
    nbin2 = morpho_count(nx,image2);
  }

  for (i=0; i<NRED; i++) imagout[i] = image1[i] & image2[i];

  if (verbose)
    message("Intersection : %d and %d -> %d\n",
            nbin1,nbin2,morpho_count(nx,imagout));

  return;
}

/*****************************************************************************/
/*!
**  Performs the union of two binary images
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  image1  first input image
** \param[in]  image2  second input image
** \param[in]  verbose Verbosity flag
**
** \param[out] imagout output image
**
*****************************************************************************/
void morpho_union(const VectorInt& nx,
                  const VectorUChar& image1,
                  const VectorUChar& image2,
                  VectorUChar& imagout,
                  bool verbose)
{
  int i,nbin1,nbin2;

  nbin1 = nbin2 = 0;
  _st_morpho_image_size_define(nx);

  if (verbose)
  {
    nbin1 = morpho_count(nx,image1);
    nbin2 = morpho_count(nx,image2);
  }

  for (i=0; i<NRED; i++) imagout[i] = image1[i] | image2[i];

  if (verbose)
    message("Union: %d and %d -> %d\n",
            nbin1,nbin2,morpho_count(nx,imagout));

  return;
}

/*****************************************************************************/
/*!
**  Performs the negation of a binary image
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  imagin  input image
** \param[in]  verbose Verbose flag
**
** \param[out] imagout output image
**
*****************************************************************************/
void morpho_negation(const VectorInt& nx,
                     const VectorUChar& imagin,
                     VectorUChar& imagout,
                     bool verbose)
{
  int i,nbin;

  nbin = 0;
  _st_morpho_image_size_define(nx);

  if (verbose) nbin = morpho_count(nx,imagin);

  for (i=0; i<NRED; i++) imagout[i] = ~imagin[i];

  if (verbose)
    message("Negation: %d -> %d\n",nbin,morpho_count(nx,imagout));

  return;
}

/*****************************************************************************/
/*!
**  Returns the volume of the grain
**
** \return  Number of pixels
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  imagin  input image
**
*****************************************************************************/
int morpho_count(const VectorInt& nx, const VectorUChar& imagin)
{
  int ix,iy,iz,ncount;

  _st_morpho_image_size_define(nx);

  ncount = 0;
  for (iz=0; iz<NX[2]; iz++)
    for (iy=0; iy<NX[1]; iy++)
      for (ix=0; ix<NX[0]; ix++)
        if (VALIN(ix,iy,iz)) ncount++;

  return(ncount);
}

/*****************************************************************************/
/*!
**  Performs a morphological opening
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  option  Option of the structuring element (CROSS pr BLOCK)
** \param[in]  radius  Radius of the structuring element (dimension = 3)
** \param[in]  imagin  input image
**
** \param[out] imagout output image
**
*****************************************************************************/
void morpho_opening(const VectorInt& nx,
                    int option,
                    const VectorInt& radius,
                    const VectorUChar& imagin,
                    VectorUChar& imagout)
{
  VectorUChar imagtmp = morpho_image_manage(nx);

  morpho_erosion (nx,option,radius,imagin,imagtmp);

  morpho_dilation(nx,option,radius,imagtmp,imagout);

  return;
}

/*****************************************************************************/
/*!
**  Performs a morphological closing
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  option  Option of the structuring element (CROSS pr BLOCK)
** \param[in]  radius  Radius of the structuring element (dimension = 3)
** \param[in]  imagin  input image
**
** \param[out] imagout output image
**
*****************************************************************************/
void morpho_closing(const VectorInt& nx,
                    int option,
                    const VectorInt& radius,
                    const VectorUChar& imagin,
                    VectorUChar& imagout)
{
  VectorUChar imagtmp = morpho_image_manage(nx);

  morpho_dilation(nx,option,radius,imagin,imagtmp);

  morpho_erosion (nx,option,radius,imagtmp,imagout);

  return;
}

/*****************************************************************************/
/*!
**  Converts an input image (double) into an image
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  tab     input array (double)
** \param[in]  vmin    Minimum value (Inclusive bound or TEST)
** \param[in]  vmax    Maximum value (Exclusive bound or TEST)
** \param[in]  verbose Verbose flag
**
** \param[out] imagout output image
**
*****************************************************************************/
void morpho_double2image(const VectorInt& nx,
                         const VectorDouble& tab,
                         double vmin,
                         double vmax,
                         VectorUChar& imagout,
                         bool verbose)
{
  int    ix,iy,iz,ind;
  unsigned char mot;
  double val,result;

  _st_morpho_image_size_define(nx);
  imagout.clear();
  mot = 0;
  ind = 0;
  for (iz=0; iz<NX[2]; iz++)
    for (iy=0; iy<NX[1]; iy++)
      for (ix=0; ix<NX[0]; ix++)
      {
        ind++;
        val = TAB(ix,iy,iz);
        result = 1.;
        if (FFFF(val)) result = 0.;
        if (! FFFF(vmin) && val <  vmin) result = 0.;
        if (! FFFF(vmax) && val >= vmax) result = 0.;
        mot = (mot << 1) + (unsigned char) result;
        if (ind == 8)
        {
          imagout.push_back(mot);
          mot = 0;
          ind = 0;
        }
      }

  if (ind != 0) imagout.push_back(mot << (8-ind));

  if (verbose)
    message("Translation: %d  / %d\n",morpho_count(nx,imagout),NXYZ);

  return;
}

/*****************************************************************************/
/*!
**  Converts an input image (double) into a returned image
**
** \returns The newly created image
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  tab     input array (double)
** \param[in]  vmin    Minimum value (Inclusive bound or TEST)
** \param[in]  vmax    Maximum value (Exclusive bound or TEST)
** \param[in]  verbose Verbose flag
**
*****************************************************************************/
VectorUChar morpho_double2image(const VectorInt& nx,
                                const VectorDouble& tab,
                                double vmin,
                                double vmax,
                                bool verbose)
{
  VectorUChar imagout = morpho_image_manage(nx);
  morpho_double2image(nx, tab, vmin, vmax, imagout, verbose);
  return imagout;
}

/*****************************************************************************/
/*!
**  Converts an image into an array (double)
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  imagin  input image
** \param[in]  mode    grain/pore assignment
** \li                  0 : if the value is assigned to tab()
** \li                  1 : if the value is added to tab()
** \li                 -1 : if the value is subtracted from tab()
** \param[in]  grain   value for the grain (double)
** \param[in]  pore    value for the pore (double)
** \param[in]  verbose Verbose flag
**
** \param[out] tab    output array (double)
**
*****************************************************************************/
void morpho_image2double(const VectorInt& nx,
                         const VectorUChar& imagin,
                         int mode,
                         double grain,
                         double pore,
                         VectorDouble& tab,
                         bool verbose)
{
  int ix,iy,iz;
  double value;

  _st_morpho_image_size_define(nx);

  if (verbose)
    message("Translation: %d / %d\n",morpho_count(nx,imagin),NXYZ);

  for (iz=0; iz<NX[2]; iz++)
    for (iy=0; iy<NX[1]; iy++)
      for (ix=0; ix<NX[0]; ix++)
      {
        value = VALIN(ix,iy,iz) ? grain : pore;
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

/*****************************************************************************/
/*!
**  Compute the distance for the grain to the edge
**
** \param[in]  nx         Number of grid meshes (dimension = 3)
** \param[in]  option     connectivity option (CROSS or BLOCK)
** \param[in]  radius     Radius of the structuring element (dimension = 3)
** \param[in]  flag_erode 1 Inflate the grain; 0 Reduce the grain
** \param[in]  imagin     input image
**
** \param[out] dist       output array containing the distances
**
*****************************************************************************/
void morpho_distance(const VectorInt& nx,
                     int option,
                     const VectorInt& radius,
                     int flag_erode,
                     VectorUChar& imagin,
                     VectorDouble& dist)
{
  int i,incr;

  /* Allocate a temporary image */

  VectorUChar imagout = morpho_image_manage(nx);

  /* Copy the initial image in the distance array */

  morpho_image2double(nx,imagin,0,1,0,dist);

  /* Processing loop */

  if (flag_erode)
  {
    while (morpho_count(nx,imagin) != 0)
    {
      morpho_erosion(nx,option,radius,imagin,imagout);
      morpho_duplicate(nx,imagout,imagin);
      morpho_image2double(nx,imagin,1,1,0,dist);
    }
  }
  else
  {
    incr = 0;
    while (morpho_count(nx,imagin) != NXYZ)
    {
      incr++;
      morpho_dilation(nx,option,radius,imagin,imagout);
      morpho_duplicate(nx,imagout,imagin);
      morpho_image2double(nx,imagin,-1,1,0,dist);
    }
    morpho_image2double(nx,imagin,1,incr,0,dist);
  }

  /* Turn the distance to a positive value */

  for (i=0; i<NXYZ; i++) dist[i] = ABS(dist[i]);

  return;
}

/*****************************************************************************/
/*!
**  Print an image stored as a bitmap
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  imagin  input image
**
*****************************************************************************/
void bitmap_print(const VectorInt& nx, const VectorUChar& imagin)
{
  int ix,iy,iz;
  unsigned char val;

  /* Initializations */

  _st_morpho_image_size_define(nx);

  /* Loop on the levels */

  for (iz=0; iz<nx[2]; iz++)
  {
    if (nx[2] > 1)
      mestitle(2,"Level %d/%d",iz+1,nx[2]);
    else
      message("\n");

    /* Loop on the cells of the layer */

    message("     ");
    for (ix=0; ix<nx[0]; ix++)
    {
      val = (ix+1) % 10;
      message("%d",val);
    }
    message("\n\n");

    for (ix=0; ix<nx[0]; ix++)
    {
      message(" %3d ",ix+1);
      for (iy=0; iy<nx[1]; iy++)
      {
        val = (VALIN(ix,iy,iz) > 0) ? 1 : 0;
        message("%d",val);
      }
      message("\n");
    }
  }
}

/*****************************************************************************/
/*!
**  Returns the size of an IMAGE
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
**
*****************************************************************************/
int bitmap_size(const VectorInt& nx)
{
  _st_morpho_image_size_define(nx);
  return(NRED);
}

/*****************************************************************************/
/*!
**  Returns the value of a bit of a bitmap array
**
** \return The value (0 or 1) of the target bit
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  imagin  Target IMAGE
** \param[in]  ix      Index of the bit along X
** \param[in]  iy      Index of the bit along Y
** \param[in]  iz      Index of the bit along Z
**
*****************************************************************************/
int bitmap_get_value(const VectorInt& nx,
                     const VectorUChar& imagin,
                     int ix,
                     int iy,
                     int iz)
{
  int retval;
  retval = (VALIN(ix,iy,iz) > 0) ? 1 : 0;
  return(retval);
}

/*****************************************************************************/
/*!
**  Set the value of a bit of a bitmap array
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  ix      Index of the bit along X
** \param[in]  iy      Index of the bit along Y
** \param[in]  iz      Index of the bit along Z
** \param[in]  imagout Target IMAGE
** \param[in]  bitval  Value of the bit (0 or 1)
**
*****************************************************************************/
void bitmap_set_value(const VectorInt& nx,
                      VectorUChar& imagout,
                      int ix,
                      int iy,
                      int iz,
                      int bitval)
{
  if (bitval > 0)
    IMAGOUT(ix,iy,iz) |= OFFSET(ix,iy,iz);
  else
    IMAGOUT(ix,iy,iz) &= MASKOFF(ix,iy,iz);
  return;
}

/****************************************************************************/
/*!
**  Create the array of index shifts for a dilation by 'radius' of a
**  regular grid
**
** \return Pointer to the newly allocated vector
**
** \param[in]  ndim        Space dimension for the grid
** \param[in]  option      Connectivity option (CROSS or BLOCK)
** \param[in]  radius      Dilation radius (0: no dilation)
** \param[in]  flag_center 1 to omit the center
** \param[in]  verbose     Verbose flag
**
** \param[out] nvois       Number of neighboring cells
**
** \remarks  The resulting array has dimension: nvois * ndim
**
*****************************************************************************/
VectorInt gridcell_neigh(int ndim,
                         int option,
                         int radius,
                         int flag_center,
                         int verbose,
                         int *nvois)
{
  int    *indg0,*indg1,ecr,flag_count,nech;
  Db     *grid;
  VectorInt    nx, indret;
  VectorDouble x0;
  VectorDouble dx;

  /* Initializations */

  (*nvois) = 0;
  indg0 = indg1 = nullptr;

  /* Create the grid attributes */

  nx.resize(ndim);
  x0.resize(ndim);
  dx.resize(ndim);
  for (int idim=0; idim<ndim; idim++)
  {
    dx[idim] = 1.;
    nx[idim] = 1 + 2 * radius;
    x0[idim] = 0.;
  }

  grid = db_create_grid(0,ndim,0,ELoadBy::SAMPLE,1,nx,x0,dx);

  /* (Maximum) core allocation */

  nech = grid->getSampleNumber();
  indret.resize(nech * ndim);
  indg0 = db_indg_alloc(grid);
  if (indg0 == nullptr) goto label_end;
  indg1 = db_indg_alloc(grid);
  if (indg1 == nullptr) goto label_end;

  /* Scan the grid nodes */

  ecr = 0;
  db_index_sample_to_grid(grid,nech/2,indg0);
  for (int iech=0; iech<nech; iech++)
  {
    db_index_sample_to_grid(grid,iech,indg1);
    flag_count = 0;
    for (int idim=0; idim<ndim; idim++)
    {
      indg1[idim] -= indg0[idim];
      if (indg1[idim] != 0) flag_count++;
    }

    if (flag_center && iech == nech/2) continue;
    if (option == CROSS && flag_count > 1) continue;

    /* Add this cell */

    for (int idim=0; idim<ndim; idim++)
      indret[ecr++] = indg1[idim];
  }

  /* Resizing the returned array */

  indret.resize(ecr);
  (*nvois) = ecr / ndim;

  /* Optional printout */

  if (verbose && (*nvois) > 0)
  {
    ecr = 0;
    message("Grid Dilation: %d samples\n",(*nvois));
    for (int i=0; i<(*nvois); i++)
    {
      message("  Neigh %3d:",i+1);
      for (int idim=0; idim<ndim; idim++)
        message(" %2d",indret[ecr++]);
      message("\n");
    }
  }

label_end:
  grid  = db_delete(grid);
  indg0 = db_indg_free(indg0);
  indg1 = db_indg_free(indg1);
  return(indret);
}

/*****************************************************************************/
/*!
**  Calculate the gradient orientations of a colored image
**
** \param[in]  nx      Number of grid meshes (dimension = 3)
** \param[in]  radius  Neighborhood dimension
** \param[in]  tab     input array (double)
**
** \param[out] tabout  Output angle
**
*****************************************************************************/
void morpho_angle(const VectorInt& nx, int radius, double *tab, double *tabout)
{
  int    ix,iy,iz,iiz,pivot,lec;
  double xi,yi,zi,z0,a[3],b[2],x[2],result;

  /* Initializations */

  _st_morpho_image_size_define(nx);
  if (NX[2] > 1)
  {
    messerr("The function is programmed for the 2-D case only");
    return;
  }

  /* Processing */

  iz = iiz = lec = 0;
  for (int iiy=0; iiy<NX[1]; iiy++)
    for (int iix=0; iix<NX[0]; iix++)
    {
      TABOUT(iix,iiy,iiz) = TEST;

      /* Collecting the neighborhood */

      for (int i=0; i<3; i++) a[i] = 0.;
      for (int i=0; i<2; i++) b[i] = x[i] = 0.;
      z0 = TAB(iix,iiy,iiz);
      if (FFFF(z0)) continue;

      for (int jx=-radius; jx<=radius; jx++)
        for (int jy=-radius; jy<=radius; jy++)
        {
          ix  = iix + jx;
          if (ix < 0 || ix >= NX[0]) continue;
          iy  = iiy + jy;
          if (iy < 0 || iy >= NX[1]) continue;
          xi  = jx;
          yi  = jy;
          zi  = TAB(ix,iy,iz);
          if (FFFF(zi)) continue;
          zi   -= z0;
          a[0] += xi * xi;
          a[1] += xi * yi;
          a[2] += yi * yi;
          b[0] += xi * zi;
          b[1] += yi * zi;
        }
      if (matrix_solve(0,a,b,x,2,1,&pivot)) continue;

      result  = ut_rad2deg(atan2(x[1],x[0]));
      result += 90.;
      while(result < -180) result += 360.;
      while(result >  180) result -= 360.;
      TABOUT(iix,iiy,iiz) = result;
    }

  return;
}

/*****************************************************************************/
/*!
**  Manage core for an image
**
** \return  Newly allocated image
**
** \param[in]  nx     Number of grid meshes (dimension = 3)
**
*****************************************************************************/
VectorUChar morpho_image_manage(const VectorInt& nx)
{
  _st_morpho_image_size_define(nx);

  VectorUChar imagout = VectorUChar(NRED);

  for (int i=0; i<NRED; i++) imagout[i] = 0;

  return (imagout);
}

