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
#include "Basic/String.hpp"
#include "Basic/Utilities.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "geoslib_old_f.h"
#include <cstring>

namespace gstlrn
{
struct SPIMG
{
  VectorDouble bitmap; /* pointer to pixels */
};

/*! \cond */

#define BORD 2

#define SURFACE_UNKNOWN 0
#define SURFACE_OUTSIDE 1
#define SURFACE_INSIDE  2
#define SURFACE_BELOW   3
#define INQUEUE         -1

#define IAD(ix, iy)          ((ix) + (iy) * TX)
#define BITMAP(im, ix, iy)   im->bitmap[(ix) + BORD + ((iy) + BORD) * TX]
#define BITALL(im, ix, iy)   im->bitmap[IAD(ix, iy)]
#define MARK(ix, iy)         *(pt_mark + IAD(ix, iy))
#define NBGH(ix, iy)         *(pt_out + IAD(ix, iy))
#define OUT_TO_IN(pt_out)    *(pt_out + Offset_out_in)
#define MARK_TO_OUT(pt_mark) *(pt_mark + Offset_mark_out)
#define OUT_TO_MARK(pt_out)  *(pt_out - Offset_mark_out)
/*! \endcond */

static Id SX; /* Window size along X */
static Id SY; /* Window size along Y */
static Id TX; /* Allocated size of BITMAP along X */
static Id TY; /* Allocated size of BITMAP along Y */
static Id TXY;
static Id SXY;
static Id STEP  = 0;
static Id Hsize = 0;

static double HMAX, HINIT, HTOP, BIGVAL;
static Id Offset_mark_out, Offset_out_in, SIGNE, OPTION, VERBOSE_STEP;
static DbGrid* DB;
static std::vector<double*> Heap;

static SPIMG* SPIMG_OUT  = nullptr;
static SPIMG* SPIMG_IN   = nullptr;
static SPIMG* SPIMG_MARK = nullptr;
static double* PT_SPILL  = nullptr;

/*****************************************************************************/
/*!
 **  Evaluates the highest elevation within the reservoir
 **
 *****************************************************************************/
static double st_htop_evaluate()
{
  double value;
  double high = -SIGNE * MAXIMUM_BIG;
  for (Id iy = 0; iy < SY; iy++)
    for (Id ix = 0; ix < SX; ix++)
    {
      if (BITMAP(SPIMG_MARK, ix, iy) != SURFACE_INSIDE) continue;
      value = BITMAP(SPIMG_IN, ix, iy);
      if (SIGNE * (value - high) > 0) high = value;
    }
  return high;
}

/*****************************************************************************/
/*!
 **  Returns the coordinates of a point, given its pointer in 'out'
 **
 ** \param[in]  pt_out      Address in the image
 ** \param[in]  image       IMAGE structure
 ** \param[in]  flag_center When TRUE, coordinates are epressed in central image
 **
 ** \param[out] ix          Location of the spill point grid node along X
 ** \param[out] iy          Location of the spill point grid node along Y

 **
 *****************************************************************************/
static void st_get_coordinates(const double* pt_out,
                               Id* ix,
                               Id* iy,
                               SPIMG* image     = SPIMG_OUT,
                               bool flag_center = false)
{
  Id shift = static_cast<Id>(pt_out - image->bitmap.data());
  *iy      = shift / TX;
  *ix      = shift % TX;

  if (flag_center)
  {
    *ix -= BORD;
    *iy -= BORD;
  }
}

/*****************************************************************************/
/*!
 **  Prints the current output flag array
 **
 ** \param[in]  flagMain TRUE if it is called from a main level
 ** \param[in]  title    Title for the dump (main level)
 ** \param[in]  pt_out   Designation of the target node in 'out' (if provided)
 ** \param[in]  image    Image containing the information to be displayed
 **
 *****************************************************************************/
static void st_dump(bool flagMain, const String& title, double* pt_out, SPIMG* image)
{
  String STRING;

  if (VERBOSE_STEP < 0) return;
  if (!flagMain && STEP <= VERBOSE_STEP) return;

  /* Process the title */

  (void)gslStrcpy(STRING, "\n");
  if (flagMain)
  {
    (void)gslSPrintfCat(STRING, "\nEnd of Step %d === ", STEP);
    (void)gslStrcat(STRING, title.c_str());
    (void)gslStrcat(STRING, "\n");
  }

  /* Current address */

  Id ix0 = ITEST;
  Id iy0 = ITEST;
  if (pt_out != nullptr)
  {
    st_get_coordinates(pt_out, &ix0, &iy0);
    (void)gslSPrintfCat(STRING, "Step %d : Node (%d/%d, %d/%d)\n", STEP, ix0, TX, iy0, TY);
  }
  message(STRING.data());

  // Current Spill position

  Id ix_spill = ITEST;
  Id iy_spill = ITEST;
  if (PT_SPILL != nullptr)
  {
    st_get_coordinates(PT_SPILL, &ix_spill, &iy_spill);
  }

  // Loop on the rows (reversed order) followed by the loop on columns
  // for a legible printout

  Id numm1 = 0;
  Id nump0 = 0;
  Id nump1 = 0;
  Id nump2 = 0;
  Id numpb = 0;
  for (Id jy = 0; jy < TY; jy++)
  {
    Id iy = TY - jy - 1;
    (void)gslStrcpy(STRING, "");
    for (Id ix = 0; ix < TX; ix++)
    {
      Id value = BITALL(image, ix, iy);
      if (ix == ix0 && iy == iy0)
      {
        (void)gslStrcat(STRING, "X");
      }
      else if (ix == ix_spill && iy == iy_spill)
      {
        (void)gslStrcat(STRING, "#");
      }
      else if (value == INQUEUE)
      {
        numm1++;
        (void)gslStrcat(STRING, "?");
      }
      else if (value == SURFACE_UNKNOWN)
      {
        nump0++;
        (void)gslStrcat(STRING, " ");
      }
      else if (value == SURFACE_OUTSIDE)
      {
        (void)gslStrcat(STRING, ".");
        nump1++;
      }
      else if (value == SURFACE_INSIDE)
      {
        (void)gslStrcat(STRING, "*");
        nump2++;
      }
      else if (value == SURFACE_BELOW)
      {
        (void)gslStrcat(STRING, "-");
        numpb++;
      }
      else
      {
        (void)gslStrcat(STRING, "U");
      }
    }
    (void)gslStrcat(STRING, "\n");
    message(STRING.data());
  }
  message("Spill'#' Queue'?'(%d) Unknown' '(%d) Out'.'(%d) In'*'(%d) Below'-'(%d) Heap(%d)\n",
          numm1, nump0, nump1, nump2, numpb, Hsize);
}

/*****************************************************************************/
/*!
 **  Blanks the center of the image
 **
 ** \param[in,out]  image SPIMG structure to be initialized
 **
 *****************************************************************************/
static void st_blank_center(SPIMG* image)
{
  for (Id iy = 0; iy < SY; iy++)
    for (Id ix = 0; ix < SX; ix++)
      BITMAP(image, ix, iy) = 0.;
}

/*****************************************************************************/
/*!
 **  Loads the center of the image from an input array
 **
 ** \param[in]  mode   Type of information
 ** \li                 0 : for the height variable
 ** \li                 1 : for the data variable
 ** \param[in]  iatt   Rank of the attribute
 ** \param[in]  defval Default value
 ** \param[in,out] image SPIMG structure to be initialized
 **
 *****************************************************************************/
static void st_copy_center(Id mode, Id iatt, SPIMG* image, double defval)
{
  DECLARE_UNUSED(mode);

  VectorInt ind(2);
  for (Id iy = 0; iy < SY; iy++)
    for (Id ix = 0; ix < SX; ix++)
    {
      ind[0]                = ix;
      ind[1]                = iy;
      double value          = DB->getArray(DB->indiceToRank(ind), iatt);
      BITMAP(image, ix, iy) = (FFFF(defval)) ? defval : value;
    }
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
static void st_extract_center(SPIMG* image, Id iatt)
{
  VectorInt ind(2);
  for (Id iy = 0; iy < SY; iy++)
    for (Id ix = 0; ix < SX; ix++)
    {
      ind[0] = ix;
      ind[1] = iy;
      DB->setArray(DB->indiceToRank(ind), iatt, BITMAP(image, ix, iy));
    }
}

static void st_change(double* pt_out, double value)
{
  *pt_out = value;
  if (value != SURFACE_UNKNOWN) st_dump(false, String(), pt_out, SPIMG_OUT);
}

/*****************************************************************************/
/*!
 **  Converts the final image into the following codes:
 **  SURFACE_INSIDE, SURFACE_OUTSIDE, SURFACE_BELOW or SURFACE_UNKNOWN
 **
 ** \param[in]  hspill  spill elevation
 **
 *****************************************************************************/
static void st_convert(double hspill)
{
  double *pt_out, th;

  for (Id iy = 0; iy < SY; iy++)
    for (Id ix = 0; ix < SX; ix++)
    {
      pt_out = &BITMAP(SPIMG_OUT, ix, iy);
      if (*pt_out == SURFACE_INSIDE || *pt_out == SURFACE_OUTSIDE) continue;

      th = SIGNE * (OUT_TO_IN(pt_out) - hspill);
      if (th < 0.)
        st_change(pt_out, SURFACE_BELOW);
      else
        st_change(pt_out, SURFACE_UNKNOWN);
    }
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
static SPIMG* st_image_free(SPIMG* image)
{
  if (image == nullptr) return (image);
  delete image;
  image = nullptr;
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
  SPIMG* image;
  Id error;

  /* Initializations */

  image = NULL;
  error = 1;

  /* Create the header */

  image = new SPIMG;

  /* Create the pixel array */

  image->bitmap.resize(TXY, value);

  /* Set the array to zero */

  error = 0;

  if (error)
  {
    delete image;
    image = nullptr;
  }
  return (image);
}

/*****************************************************************************/
/*!
 **  Add an element to the Heap Sort Pile
 **
 ** \param[in]  p pointer to the element to be added
 **
 *****************************************************************************/
static void st_heap_add(double* p)
{
  Id i    = Hsize++;
  Id n    = (i - 1) / 2;
  Heap[i] = p;
  while ((i > 0) && SIGNE * (OUT_TO_IN(p) - OUT_TO_IN(Heap[n])) > 0.)
  {
    Heap[i] = Heap[n];
    i       = n;
    n       = (i - 1) / 2;
  }

  Heap[i] = p;
  *p      = INQUEUE;
  st_dump(false, String(), p, SPIMG_OUT);
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
  Id i, il, ir, is;

  first = Heap[0];
  Hsize--;
  Heap[0] = Heap[Hsize];
  i       = 0;
  while (i < Hsize / 2)
  {
    il = 2 * i + 1;
    ir = 2 * i + 2;
    is = i;
    if ((il < Hsize) && SIGNE * (OUT_TO_IN(Heap[il]) - OUT_TO_IN(Heap[is])) > 0)
      is = il;
    if ((ir < Hsize) && SIGNE * (OUT_TO_IN(Heap[ir]) - OUT_TO_IN(Heap[is])) > 0)
      is = ir;

    if (is == i) break;
    temp     = Heap[i];
    Heap[i]  = Heap[is];
    Heap[is] = temp;
    i        = is;
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
static Id st_traite(double* pt_out, double* pt_vois)

{
  double th, value;

  if (*pt_vois == SURFACE_OUTSIDE || *pt_vois == SURFACE_INSIDE)
  {

    /* The neighboring element is defined */

    if (*pt_out == INQUEUE)
    {

      /* Copy the value of the neighboring element */

      st_change(pt_out, *pt_vois);
      if (*pt_vois == SURFACE_INSIDE)
      {
        value = OUT_TO_IN(pt_out);
        if (SIGNE * (value - HTOP) > 0) HTOP = value;
        if (!FFFF(HMAX))
        {
          th = SIGNE * (HTOP - value);
          if (th > HMAX) return (2);
        }
      }
    }
    else if (*pt_out != SURFACE_UNKNOWN && *pt_out != *pt_vois)
    {

      /* The neighbor is already valued but differently */

      st_change(pt_out, SURFACE_INSIDE);
      value = OUT_TO_IN(pt_out);
      if (SIGNE * (value - HTOP) > 0) HTOP = value;
      if (!FFFF(HMAX))
      {
        th = SIGNE * (HTOP - value);
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

static void st_print()
{
  mestitle(1, "Spill Point environment");
  message("- Grid dimensions = %d x %d\n", SX, SY);
  if (!FFFF(HMAX))
    message("- Maximum reservoir thickness = %lf\n", HMAX);
  else
    message("- No Maximum reservoir thickness\n");
  if (OPTION == 0)
    message("- 4 - connectivity\n");
  else
    message("- 8 - connectivity\n");
  message("An edge of %d pixels is added to the original grid.\n", BORD);
}

static void st_final_stats(double hspill, Id ix0, Id iy0)
{
  Id num_inside      = 0;
  double min_inside  = MAXIMUM_BIG;
  double max_inside  = MINIMUM_BIG;
  Id num_outside     = 0;
  double min_outside = MAXIMUM_BIG;
  double max_outside = MINIMUM_BIG;
  Id num_else        = 0;
  double min_else    = MAXIMUM_BIG;
  double max_else    = MINIMUM_BIG;

  for (Id iy = 0; iy < SY; iy++)
    for (Id ix = 0; ix < SX; ix++)
    {
      Id value    = BITMAP(SPIMG_OUT, ix, iy);
      double topo = BITMAP(SPIMG_IN, ix, iy);

      if (value == SURFACE_INSIDE)
      {
        if (topo < min_inside) min_inside = topo;
        if (topo > max_inside) max_inside = topo;
        num_inside++;
      }
      else if (value == SURFACE_OUTSIDE)
      {
        if (topo < min_outside) min_outside = topo;
        if (topo > max_outside) max_outside = topo;
        num_outside++;
      }
      else
      {
        if (topo < min_else) min_else = topo;
        if (topo > max_else) max_else = topo;
        num_else++;
      }
    }

  mestitle(1, "Final statistics");
  message("INSIDE:  Topography within [%lf ; %lf] (%d)\n",
          min_inside, max_inside, num_inside);
  message("OUTSIDE: Topography within [%lf ; %lf] (%d)\n",
          min_outside, max_outside, num_outside);
  message("UNKNOWN: Topography within [%lf ; %lf] (%d)\n",
          min_else, max_else, num_else);
  message("Elevation: HINIT = %lf - Spill = %lf\n", HINIT, hspill);
  message("Grid indices of the Spill Point = %d %d\n", ix0, iy0);
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
 ** \param[in]  option        0 for 4-connectivity; 1 for 8-connectivity
 ** \param[in]  flag_up       TRUE when working in elevation; 0 in depth
 ** \param[in]  verbose_step  Step for verbose flag
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
Id spill_point(DbGrid* dbgrid,
               Id ind_depth,
               Id ind_data,
               Id option,
               bool flag_up,
               Id verbose_step,
               double hmax,
               double* h,
               const double* th,
               Id* ix0,
               Id* iy0)
{
  DECLARE_UNUSED(th);
  double *pt_mark, *pt_out, hspill;
  Id *x, *y, k, n, iy, ix, found, local;
  static Id n4   = 4;
  static Id n8   = 8;
  static Id x4[] = {1, -1, 0, 0};
  static Id y4[] = {0, 0, 1, -1};
  static Id x8[] = {1, -1, 0, 0, 1, -1, -1, 1};
  static Id y8[] = {0, 0, 1, -1, 1, -1, 1, -1};

  /* Preliminary tests */

  Id error = 1;

  /* Preliminary checks */

  if (!dbgrid->isGrid())
  {
    messerr("The Spill Point algorithm is restricted to regular grid");
    return (1);
  }
  if (dbgrid->getNDim() != 2)
  {
    messerr("Spill point is limited to 2-D space");
    return (1);
  }
  if (ind_depth < 0 || ind_depth > dbgrid->getNColumn() ||
      ind_data < 0 || ind_data > dbgrid->getNColumn())
  {
    messerr("Error in the ranks of the height (%d) and data (%d) variables",
            ind_depth, ind_data);
    return (1);
  }

  /* Define global variables */

  hspill       = TEST;
  HMAX         = hmax;
  SIGNE        = (flag_up) ? 1 : -1;
  BIGVAL       = (flag_up) ? MAXIMUM_BIG : MINIMUM_BIG;
  VERBOSE_STEP = verbose_step;
  OPTION       = option;
  DB           = dbgrid;
  SX           = DB->getNX(0);
  SY           = DB->getNX(1);
  TX           = SX + 2 * BORD;
  TY           = SY + 2 * BORD;
  TXY          = TX * TY;
  SXY          = SX * SY;

  /* Initializations */

  if (OPTION)
  {
    n = n8;
    x = x8;
    y = y8;
  }
  else
  {
    n = n4;
    x = x4;
    y = y4;
  }
  pt_out = pt_mark = NULL;

  // Print the environment

  st_print();

  /* Add the attribute */

  Id iptr_spill = dbgrid->addColumnsByConstant(1, 0., "Spill", ELoc::Z);
  if (iptr_spill < 0) goto label_end;

  /* Core allocation */

  SPIMG_IN = st_image_alloc(BIGVAL);
  if (SPIMG_IN == nullptr) goto label_end;
  SPIMG_MARK = st_image_alloc(SURFACE_OUTSIDE);
  if (SPIMG_MARK == nullptr) goto label_end;
  SPIMG_OUT = st_image_alloc(SURFACE_OUTSIDE);
  if (SPIMG_OUT == nullptr) goto label_end;
  Offset_out_in   = static_cast<Id>(SPIMG_IN->bitmap.data() - SPIMG_OUT->bitmap.data());
  Offset_mark_out = static_cast<Id>(SPIMG_OUT->bitmap.data() - SPIMG_MARK->bitmap.data());

  /* Copying the input arrays into the corresponding images */

  STEP = 1;
  st_copy_center(0, ind_depth, SPIMG_IN, BIGVAL);
  st_copy_center(1, ind_data, SPIMG_MARK, SURFACE_OUTSIDE);
  st_blank_center(SPIMG_OUT);
  st_dump(true, "Constraints", NULL, SPIMG_MARK);

  HTOP = HINIT = st_htop_evaluate();

  /* Creation of the Heap-search Pile */

  Hsize = 0;
  Heap.resize(TXY);

  /***************************/
  /* Add markers to the Heap */
  /***************************/

  STEP = 2;
  for (iy = 0; iy < SY; iy++)
    for (ix = 0; ix < SX; ix++)
    {
      pt_mark = &BITMAP(SPIMG_MARK, ix, iy);
      pt_out  = &MARK_TO_OUT(pt_mark);
      if (static_cast<Id>(MARK(0, 0)) == SURFACE_UNKNOWN)
      {
        for (k = found = 0; k < n && found == 0; k++)
          if (static_cast<Id>(MARK(x[k], y[k])) == SURFACE_INSIDE) found = 1;
        if (found)
          st_heap_add(pt_out);
        else
          st_change(pt_out, SURFACE_UNKNOWN);
      }
      else if (static_cast<Id>(MARK(0, 0)) == SURFACE_INSIDE)
        st_change(pt_out, SURFACE_INSIDE);
      else
        st_heap_add(pt_out);
    }
  st_dump(true, "Markers posted", NULL, SPIMG_OUT);

  /***************/
  /* Propagation */
  /***************/

  STEP  = 3;
  found = 0;
  while (Hsize > 0)
  {
    pt_out  = st_heap_del();
    pt_mark = &OUT_TO_MARK(pt_out);
    if (*pt_mark == SURFACE_OUTSIDE) *pt_out = SURFACE_OUTSIDE;
    for (k = 0; k < n && found == 0; k++)
    {
      local = st_traite(pt_out, &NBGH(x[k], y[k]));
      found = MAX(found, local);
    }
    if (found == 0) PT_SPILL = pt_out;
  }
  st_dump(true, "After Propagation", NULL, SPIMG_OUT);

  if (found == 2 || PT_SPILL == nullptr)
  {
    hspill   = (SIGNE > 0) ? HTOP - HMAX : HMAX + HTOP;
    PT_SPILL = pt_out;
  }
  else
    hspill = OUT_TO_IN(PT_SPILL);
  st_get_coordinates(PT_SPILL, ix0, iy0, SPIMG_OUT, true);

  /******************************************************/
  /* Fill the remaining part of the flat Spill boundary */
  /******************************************************/

  STEP = 4;
  while (Hsize > 0)
  {
    pt_out = st_heap_del();
    if (OUT_TO_IN(pt_out) != hspill) break;
    for (k = 0; k < n; k++)
      (void)st_traite(pt_out, &NBGH(x[k], y[k]));
    st_change(pt_out, SURFACE_INSIDE);
    st_dump(false, String(), pt_out, SPIMG_OUT);
  }
  Heap.clear();
  st_dump(true, "After Filling Flat Boundaries", NULL, SPIMG_OUT);

  /***********************************/
  /* Final conversions of the output */
  /***********************************/

  STEP = 5;
  st_convert(hspill);
  st_dump(true, "After Final Conversion", NULL, SPIMG_OUT);

  // Print final statistics

  st_final_stats(hspill, *ix0, *iy0);

  /* Returning the output grid */

  st_extract_center(SPIMG_OUT, iptr_spill);

label_end:
  SPIMG_IN   = st_image_free(SPIMG_IN);
  SPIMG_OUT  = st_image_free(SPIMG_OUT);
  SPIMG_MARK = st_image_free(SPIMG_MARK);

  /* Returning arguments */

  *h = hspill;

  return (error);
}
} // namespace gstlrn
