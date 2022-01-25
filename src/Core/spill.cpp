/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"

#include <string.h>

typedef struct
{
  double *bitmap; /* pointer to pixels */
} SPIMG;

/*! \cond */
#define BORD                        2
#define SURFACE_UNKNOWN             0
#define SURFACE_OUTSIDE             1
#define SURFACE_INSIDE              2
#define SURFACE_BELOW               3
#define INQUEUE                    -1

#define IAD(ix,iy)         ((iy) + (ix) * TY)
#define BITMAP(im,ix,iy)  *(im->bitmap + (iy) + BORD + ((ix) + BORD) * TY)
#define BITALL(im,ix,iy)  *(im->bitmap + IAD(ix,iy))
#define MARK(ix,iy)       *(pt_mark    + IAD((ix),(iy)))
#define NBGH(ix,iy)       *(pt_out     + IAD((ix),(iy)))
#define IN_OUT(pt_out)    *(pt_out     + Offset_out_in)
#define MARK_OUT(pt_mark) *(pt_mark    + Offset_mark_out)
#define OUT_MARK(pt_out)  *(pt_out     - Offset_mark_out)
/*! \endcond */

static int SX; /* Window size along X */
static int SY; /* Window size along Y */
static int TX; /* Allocated size of BITMAP along X */
static int TY; /* Allocated size of BITMAP along Y */

static double **Heap, HMAX, HTOP;
static int Hsize, Offset_mark_out, Offset_out_in;
static int SIGNE, FLAG_VERBOSE, FLAG_CROSS;
static Db *DB;

/*****************************************************************************/
/*!
 **  Evaluates the highest elevation within the reservoir
 **
 ** \param[in]  data SPIMG structure to be initialized
 **
 ** \param[out] input SPIMG structure to be initialized
 **
 *****************************************************************************/
static void st_htop_evaluate(SPIMG *data, SPIMG *input)
{
  double value;
  int ix, iy;

  HTOP = -SIGNE * 1.e30;

  for (ix = 0; ix < SX; ix++)
    for (iy = 0; iy < SY; iy++)
    {
      if (BITMAP(data,ix,iy) != SURFACE_INSIDE) continue;
      value = BITMAP(input, ix, iy);
      if (SIGNE * (value - HTOP) > 0) HTOP = value;
    }
  return;
}

/*****************************************************************************/
/*!
 **  Blanks the center of the image
 **
 ** \param[in,out]  image SPIMG structure to be initialized
 **
 *****************************************************************************/
static void st_blank_center(SPIMG *image)
{
  int ix, iy;

  for (ix = 0; ix < SX; ix++)
    for (iy = 0; iy < SY; iy++)
      BITMAP(image,ix,iy) = 0.;

  return;
}

/*****************************************************************************/
/*!
 **  Loads the center of the image from an input array
 **
 ** \param[in]  mode   Type of information
 ** \li                 0 : for the height variable
 ** \li                 1 : for the data variable
 ** \param[in]  iatt   Rank of the attribute
 ** \param[in,out] image SPIMG structure to be initialized
 **
 *****************************************************************************/
static void st_copy_center(int mode, int iatt, SPIMG *image)
{
  int ix, iy, lec, n_in, n_out;
  double value;

  n_in = n_out = 0;
  for (ix = lec = 0; ix < SX; ix++)
    for (iy = 0; iy < SY; iy++, lec++)
    {
      value = DB->getArray(lec, iatt);

      if (mode)
      {
        if (value == SURFACE_INSIDE)
          n_in++;
        else if (value == SURFACE_OUTSIDE)
          n_out++;
        else
          value = 0.;
      }
      BITMAP(image,ix,iy) = value;
    }

  /* Optional printout */

  if (mode == 1 && FLAG_VERBOSE)
  {
    message("Conditioning data:\n");
    message("- Number of nodes inside the reservoir  = %d\n", n_in);
    message("- Number of nodes outside the reservoir = %d\n", n_out);
  }

  return;
}

/*****************************************************************************/
/*!
 **  Extracts an output array from the the center of the image
 **
 ** \param[in]  image SPIMG structure
 **
 ** \param[out] iatt  Rank of the output attribute
 **
 *****************************************************************************/
static void st_extract_center(SPIMG *image, int iatt)
{
  int ix, iy, lec;

  for (ix = lec = 0; ix < SX; ix++)
    for (iy = 0; iy < SY; iy++, lec++)
      DB->setArray(lec, iatt, BITMAP(image, ix, iy));
  return;
}

/*****************************************************************************/
/*!
 **  Converts the final image into the following codes: SURFACE_BELOW,
 **  SURFACE_INSIDE, SURFACE_OUTSIDE, SURFACE_UNKNOWN
 **
 ** \param[in,out] out  SPIMG structure
 ** \param[in]  hspill  spill evelation
 **
 *****************************************************************************/
static void st_convert(SPIMG *out, double hspill)
{
  double *pt_out, th;
  int ix, iy;

  /* Convert the SPIMG */

  for (ix = 0; ix < SX; ix++)
    for (iy = 0; iy < SY; iy++)
    {
      pt_out = &BITMAP(out, ix, iy);
      th = (SIGNE > 0) ? IN_OUT(pt_out) - hspill :
                         hspill - IN_OUT(pt_out);
      if (th < 0.)
      {
        *pt_out = SURFACE_BELOW;
      }
      else
      {
        if (*pt_out == SURFACE_INSIDE || *pt_out == SURFACE_OUTSIDE) continue;
        *pt_out = SURFACE_UNKNOWN;
      }
    }
  return;
}

/*****************************************************************************/
/*!
 **  Converts the UNKNOWN into INSIDE Check if the Maximum Reservoir
 **  Thickness criterion is still honored
 **
 ** \return  Error returned code: 1 if the Maximum Reservoir Thickness
 ** \return  constraint is violated; 0 otherwise
 **
 ** \param[in,out] out SPIMG structure
 ** \param[in]  hspill Spill elevation value
 **
 *****************************************************************************/
static int st_add_unknown(SPIMG *out, double hspill)
{
  double *pt_out, th_max, th;
  int ix, iy;

  /* Convert Unknown into Inside */

  th_max = 0.;
  for (ix = 0; ix < SX; ix++)
    for (iy = 0; iy < SY; iy++)
    {
      pt_out = &BITMAP(out, ix, iy);
      if (*pt_out != SURFACE_UNKNOWN) continue;
      th = (SIGNE > 0) ? IN_OUT(pt_out) - hspill :
                         hspill - IN_OUT(pt_out);
      if (th < 0.) continue;
      *pt_out = SURFACE_INSIDE;
      if (th > th_max) th_max = th;
    }

  /* If a Maximum Reservoir Thickness is used, check if the constraints */
  /* is still honored */

  if (FFFF(HMAX) && th_max > HMAX) return (1);

  return (0);
}

/*****************************************************************************/
/*!
 **  Procedure to free an already existig image
 **
 ** \return  Pointer to the image freed (NULL)
 **
 ** \param[in]  image pointer to the image to be freed
 **
 *****************************************************************************/
static SPIMG* st_image_free(SPIMG *image)

{
  if (image == (SPIMG*) NULL) return (image);

  if (image->bitmap != nullptr)
    image->bitmap = (double*) mem_free((char* ) image->bitmap);
  image = (SPIMG*) mem_free((char* ) image);

  return (image);
}

/*****************************************************************************/
/*!
 **  Procedure to allocates a NEW image The bitmap of the new image is
 **  set to zero
 **
 ** \return  Pointer to the new image
 **
 ** \param[in]  value conventional value for initialization
 **
 *****************************************************************************/
static SPIMG* st_image_alloc(double value)
{
  SPIMG *image;
  double *pt;
  int i, error;

  /* Initializations */

  image = NULL;
  error = 1;

  /* Create the header */

  image = (SPIMG*) mem_alloc(sizeof(SPIMG), 0);
  if (image == (SPIMG*) NULL) goto label_end;

  /* Create the pixel array */

  image->bitmap = (double*) mem_alloc(sizeof(double) * TX * TY, 0);
  if (image->bitmap == nullptr) goto label_end;

  /* Set the array to zero */

  pt = image->bitmap;
  for (i = 0; i < TX * TY; i++)
    *pt++ = value;
  error = 0;

  label_end: if (error) image = st_image_free(image);
  return (image);
}

/*****************************************************************************/
/*!
 **  Add an element to the Heap Sort Pile
 **
 ** \param[in]  p pointer to the element to be added
 **
 *****************************************************************************/
static void st_heap_add(double *p)

{
  int i, n;

  i = Hsize++;
  n = (i - 1) / 2;
  Heap[i] = p;
  while ((i > 0) && SIGNE * (IN_OUT(p) - IN_OUT(Heap[n])) > 0.)
  {
    Heap[i] = Heap[n];
    i = n;
    n = (i - 1) / 2;
  }

  Heap[i] = p;
  *p = INQUEUE;

  return;
}

/*****************************************************************************/
/*!
 **  Return the first element of the Heap Sort Pile and delete it
 **
 ** \return  Pointer to the first element of the Heap Sort Pile
 **
 *****************************************************************************/
static double* st_heap_del(void)

{
  double *first, *temp;
  int i, il, ir, is;

  first = Heap[0];
  Hsize--;
  Heap[0] = Heap[Hsize];
  i = 0;
  while (i < Hsize / 2)
  {
    il = 2 * i + 1;
    ir = 2 * i + 2;
    is = i;
    if ((il < Hsize) && SIGNE * (IN_OUT(Heap[il]) - IN_OUT(Heap[is])) > 0)
      is = il;
    if ((ir < Hsize) && SIGNE * (IN_OUT(Heap[ir]) - IN_OUT(Heap[is])) > 0)
      is = ir;

    if (is == i) break;
    temp = Heap[i];
    Heap[i] = Heap[is];
    Heap[is] = temp;
    i = is;
  }

  return (first);
}

/*****************************************************************************/
/*!
 **  Checks if a current element can be processed according to its
 **  neighborhood status
 **
 ** \return  Flag indicating the end of the procedure:
 ** \return  1 : the current element and its neighboring one have
 ** \remark      two different status
 ** \remark  2 : the maximum reservoir thickness has been reached
 **
 ** \param[in]  pt_out  pointer to the current element
 ** \param[in]  pt_vois pointer to the neighboring element
 **
 *****************************************************************************/
static int st_traite(double *pt_out, double *pt_vois)

{
  double th, value;

  if (*pt_vois == SURFACE_OUTSIDE || *pt_vois == SURFACE_INSIDE)
  {

    /* The neighboring element is defined */

    if (*pt_out == INQUEUE)
    {

      /* Copy the value of the neighboring element */

      *pt_out = *pt_vois;
      if (*pt_vois == SURFACE_INSIDE)
      {
        value = IN_OUT(pt_out);
        if (SIGNE * (value - HTOP) > 0) HTOP = value;
        if (!FFFF(HMAX))
        {
          th = (SIGNE > 0) ? HTOP - value :
                             value - HTOP;
          if (th > HMAX) return (2);
        }
      }
    }
    else if (*pt_out != SURFACE_UNKNOWN && *pt_out != *pt_vois)
    {

      /* The neighbor is already valuated but differently */

      *pt_out = SURFACE_INSIDE;
      value = IN_OUT(pt_out);
      if (SIGNE * (value - HTOP) > 0) HTOP = value;
      if (!FFFF(HMAX))
      {
        th = (SIGNE > 0) ? HTOP - value :
                           value - HTOP;
        if (th > HMAX) return (2);
      }
      return (1);
    }
  }
  else if ((*pt_vois) == SURFACE_UNKNOWN)
  {

    /* The neighboring element is added to the Heap Sort Pile */

    st_heap_add(pt_vois);
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Prints the current output flag array
 **
 ** \param[in]  title   title for the conditional dump
 ** \param[in]  pt_out  Pointer to the array to be printed
 ** \param[in]  out     working image containing the marked zones
 **
 *****************************************************************************/
static void st_dump(const char *title, double *pt_out, SPIMG *out)

{
  int ix, iy, shift;
  char STRING[BUFFER_LENGTH];
  static int first = 0;

  if (!FLAG_VERBOSE) return;

  /* Print the general comment */

  if (first == 0)
  {
    message("For the DEBUGGING print, we use the following conventions\n"
            "Pixels are organized by increasing IX in a line\n"
            "                     by increasing IY in a column\n"
            "The origin (IX=0,IY=0) is located in the upper left corner\n");
    message("A frame (%d pixels) is added to the original grid\n\n", BORD);
    first = 1;
  }

  /* Process the title */

  (void) gslStrcpy(STRING, title);
  (void) gslStrcat(STRING, " : ");

  /* Process the address */

  if (pt_out != nullptr)
  {
    shift = static_cast<int>(pt_out - out->bitmap);
    ix = shift / TY;
    iy = shift - TY * ix;
    (void) gslSPrintf(&STRING[strlen(STRING)], "Node processed=(%d, %d)", ix,
                      iy);
  }

  /* Print the header (if any) */

  (void) gslStrcat(STRING, "\n");
  message(STRING);

  /* Dump the grid */

  for (iy = 0; iy < TY; iy++)
  {
    for (ix = 0; ix < TX; ix++)
      message(" %4d", BITALL(out, ix, iy));
    message("\n");
  }
  message("\n");

  return;
}

/*****************************************************************************/
/*!
 **  Returns the coordinates of a point, given its pointer
 **
 ** \param[in]  pt_out       pointer to the point of interest
 ** \param[in]  out          IMAGE structure
 **
 ** \param[out] ix0 location of the spill point grid node along X (The numbering
 **                 must start with 1)
 ** \param[out] iy0 location of the spill point grid node along Y (The numbering
 **                 must start with 1)
 **
 *****************************************************************************/
static void st_get_coordinates(double *pt_out, SPIMG *out, int *ix0, int *iy0)
{
  *ix0 = static_cast<int>((pt_out - out->bitmap) / TY - BORD + 1);
  *iy0 = static_cast<int>((pt_out - out->bitmap) % TY - BORD + 1);

  if (OptDbg::query(EDbg::MORPHO))
    message("Processed grid node : IX=%d IY=%d - Status=%d\n", (*ix0), (*iy0),
            (*pt_out));
  return;
}

/*****************************************************************************/
/*!
 **  Establishes the spill point
 **
 ** \return  Error return code
 **
 ** \param[in]  in         image containing the input variable of interest
 ** \param[in]  mark       image containing the constraining markers
 **
 ** \param[out] out working image containing the marked zones
 ** \param[out] h   elevation of the spill point
 ** \param[out] ix0 location of the spill point grid node along X (The numbering
 **                 must start with 1)
 ** \param[out] iy0 location of the spill point grid node along Y (The numbering
 **                 must start with 1)
 **
 *****************************************************************************/
static int st_spill(SPIMG *in,
                    SPIMG *mark,
                    SPIMG *out,
                    double *h,
                    int *ix0,
                    int *iy0)
{
  double *pt_mark, *pt_out, hspill;
  int *x, *y, k, n, iy, ix, found, local;
  static int n4 = 4;
  static int n8 = 8;
  static int x4[] = { 1, -1, 0, 0 };
  static int y4[] = { 0, 0, 1, -1 };
  static int x8[] = { 1, -1, 0, 0, 1, -1, -1, 1 };
  static int y8[] = { 0, 0, 1, -1, 1, -1, 1, -1 };

  /* Initializations */

  found = 0;
  pt_out = NULL;
  Offset_out_in = static_cast<int>(in->bitmap - out->bitmap);
  Offset_mark_out = static_cast<int>(out->bitmap - mark->bitmap);
  if (FLAG_CROSS)
  {
    n = n4;
    x = x4;
    y = y4;
  }
  else
  {
    n = n8;
    x = x8;
    y = y8;
  }

  /* Creation of the Heap-search Pile */

  Hsize = 0;
  Heap = (double**) mem_alloc(sizeof(double*) * TX * TY, 0);
  if (Heap == nullptr) return (1);

  /* Initialie stage */

  st_dump("Depth/Elevation Map", NULL, in);
  st_dump("Constraints", NULL, mark);

  /***************************/
  /* Add markers to the Heap */
  /***************************/

  for (ix = -1; ix <= SX; ix++)
    for (iy = -1; iy <= SY; iy++)
    {
      pt_mark = &BITMAP(mark, ix, iy);
      pt_out = &MARK_OUT(pt_mark);
      if (MARK(0,0) == SURFACE_UNKNOWN)
      {
        for (k = found = 0; k < n && found == 0; k++)
          if (MARK(x[k],y[k]) == SURFACE_INSIDE) found = 1;
        if (found)
          st_heap_add(pt_out);
        else
          pt_out = SURFACE_UNKNOWN;
      }
      else if (MARK(0,0) == SURFACE_INSIDE)
        *pt_out = MARK(0, 0);
      else
        st_heap_add(pt_out);
    }
  st_dump("Markers posted", NULL, out);

  /***************/
  /* Propagation */
  /***************/

  while (Hsize > 0)
  {
    pt_out = st_heap_del();
    pt_mark = &OUT_MARK(pt_out);
    if (*pt_mark == SURFACE_OUTSIDE) *pt_out = SURFACE_OUTSIDE;
    for (k = found = 0; k < n && found == 0; k++)
    {
      local = st_traite(pt_out, &NBGH(x[k], y[k]));
      found = MAX(found, local);
    }
    if (found) break;
    st_dump("Propagation (4-connectivity)", pt_out, out);
  }

  /* Process the last cell */

  if (found == 2)
    hspill = (SIGNE > 0) ? HTOP - HMAX :
                           HMAX + HTOP;
  else
    hspill = IN_OUT(pt_out);
  st_get_coordinates(pt_out, out, ix0, iy0);

  /******************************************************/
  /* Fill the remaining part of the flat Spill boundary */
  /******************************************************/

  while (Hsize > 0)
  {
    pt_out = st_heap_del();
    if (IN_OUT(pt_out) != hspill) break;
    for (k = 0; k < n; k++)
      (void) st_traite(pt_out, &NBGH(x[k], y[k]));
    *pt_out = SURFACE_INSIDE;
    st_dump("Filling Flat Boundary (4-connectivity)", pt_out, out);
  }

  /*****************************************/
  /* Final conversions of the output SPIMG */
  /*****************************************/

  st_convert(out, hspill);

  /* Core deallocation */

  Heap = (double**) mem_free((char* ) Heap);

  /* Returning argument */

  *h = hspill;

  return (0);
}

/*****************************************************************************/
/*!
 **  Evaluates the spill point
 **
 ** \return  Error return code
 ** \return  - Memory problem
 ** \return  - Maximum Reservoir Thickness violation when turning UNKNOWN into
 ** \return  INSIDE
 **
 ** \param[in]  dbgrid        Grid Db structure
 ** \param[in]  ind_depth     Rank of the variable containing the depth
 ** \param[in]  ind_data      Rank of the variable containing the data
 ** \param[in]  flag_up       1 when working in elevation; 0 in depth
 ** \param[in]  flag_cross    1 for 4-connectivity; 0 for 8-connectivity
 ** \param[in]  flag_unknown  1 if Unknown must be converted into Reservoir;
 **                           0 otherwise
 ** \param[in]  flag_verbose  1 for a verbose output
 ** \param[in]  hmax          maximum reservoir thickness (FFFF not used)
 **
 ** \param[out] h      elevation of the spill point
 ** \param[out] th     maximum reservoir thickness
 ** \param[out] ix0    location of the spill point grid node along X
 ** \param[out] iy0    location of the spill point grid node along Y
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
int spill_point(Db *dbgrid,
                                int ind_depth,
                                int ind_data,
                                int flag_up,
                                int flag_cross,
                                int flag_unknown,
                                int flag_verbose,
                                double hmax,
                                double *h,
                                double *th,
                                int *ix0,
                                int *iy0)
{
  SPIMG *in, *out, *mark;
  double hspill, thick;
  int error, iptr_spill;

  /* Preliminary tests */

  error = 1;

  /* Preliminary checks */

  if (!is_grid(dbgrid))
  {
    messerr("The Fluid Propagation is restricted to regular grid");
    return (1);
  }
  if (dbgrid->getNDim() != 2)
  {
    messerr("Spill point is limited to 2-D space");
    return (1);
  }
  if (ind_depth < 0 || ind_depth > dbgrid->getFieldNumber() || ind_data < 0
      || ind_data > dbgrid->getFieldNumber())
  {
    messerr("Error in the ranks of the height (%d) and data (%d) variables",
            ind_depth, ind_data);
    return (1);
  }

  /* Define global variables */

  thick = 0.;
  hspill = TEST;
  HMAX = hmax;
  SIGNE = (flag_up) ? 1 :
                      -1;
  FLAG_VERBOSE = flag_verbose;
  FLAG_CROSS = flag_cross;
  DB = dbgrid;
  SX = DB->getNX(0);
  SY = DB->getNX(1);
  TX = SX + 2 * BORD;
  TY = SY + 2 * BORD;
  in = out = mark = (SPIMG*) NULL;

  /* Add the attribute */

  iptr_spill = dbgrid->addFieldsByConstant(1, 0.);
  if (iptr_spill < 0) goto label_end;

  /* Core allocation */

  in = st_image_alloc(SURFACE_UNKNOWN);
  if (in == (SPIMG*) NULL) goto label_end;
  mark = st_image_alloc(SURFACE_OUTSIDE);
  if (mark == (SPIMG*) NULL) goto label_end;
  out = st_image_alloc(SURFACE_OUTSIDE);
  if (out == (SPIMG*) NULL) goto label_end;

  /* Copying the input arrays into the corresponding images */

  st_copy_center(0, ind_depth, in);
  st_copy_center(1, ind_data, mark);
  st_blank_center(out);
  st_htop_evaluate(mark, in);

  /* Calling the Spill calculation routine */

  error = st_spill(in, mark, out, &hspill, ix0, iy0);

  /* Convert UNKNOWN into INSIDE (upon request) */

  if (flag_unknown) error = st_add_unknown(out, hspill);
  thick = (SIGNE > 0) ? HTOP - hspill :
                        hspill - HTOP;

  /* Returning the output grid */

  st_extract_center(out, iptr_spill);

  label_end:

  /* Core deallocation */

  in = st_image_free(in);
  out = st_image_free(out);
  mark = st_image_free(mark);

  /* Returning arguments */

  *h = hspill;
  *th = thick;

  return (error);
}
