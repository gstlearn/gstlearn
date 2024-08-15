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
#include "geoslib_f_private.h"

#include "Mesh/MeshETurbo.hpp"
#include "LinearOp/ShiftOpCs.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/OptimCostColored.hpp"
#include "Stats/Classical.hpp"
#include "Covariances/CovAniso.hpp"
#include "Model/ANoStat.hpp"
#include "Model/Model.hpp"
#include "Model/CovInternal.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbHelper.hpp"
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/PolyLine2D.hpp"
#include "Polygon/Polygons.hpp"
#include "Tree/Ball.hpp"

#include <math.h>
#include <string.h>

// https://stackoverflow.com/a/26359433/3952924
#ifdef _MSC_VER
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

/*! \cond */
#define TRACE(i,iseg)       (trace[(i) * nseg + (iseg)])
#define LINE(nbline,i)      (line[npline * (nbline) + (i)])
#define PROP1(iz,iprop)     (prop1[(iz) * nprop + (iprop)])
#define PROP2(iz,iprop)     (prop2[(iz) * nprop + (iprop)])
#define WTAB(iz,icode,ivar) (wtab[(ivar) + nvar * ((iz) + nz * (icode))])
#define WCOR(iz,icode,idim) (wcor[(idim) + ndim * ((iz) + nz * (icode))])
#define WCNT(iz,icode)      (wcnt[                 (iz) + nz * (icode)])

#define R(i,j)              (R[(i) * n + (j)])

typedef struct
{
  char key[5];
  char title[STRING_LENGTH];
  int flag_rank;
  int flag_bounds;
} Edit_Item;
static int N_EDIT = 10;
static Edit_Item EDIT[] = { { "P", "Define the Properties", 0, 0 },
                            { ".", "Same as previous Command", 0, 0 },
                            { "S", "Shift the Sample Rank", 1, 0 },
                            { "V", "Shift the Variable Rank", 1, 0 },
                            { "AS", "Set the Sample Rank", 1, 0 },
                            { "AV", "Set the Variable Rank", 1, 0 },
                            { "M", "Modify the current Value", 0, 0 },
                            { "FD", "Find Next value in Interval", 0, 1 },
                            { "FU", "Find Next value in Interval", 0, 1 },
                            { "D", "Current Display", 0, 0 } };

/*! \endcond */

/****************************************************************************/
/*!
 **  Calculate the (discretized) surface of influence
 **
 ** \return  Error returned code
 **
 ** \param[in]  db_point Db containing the data points
 ** \param[in]  db_grid  Db containing the discretization grid
 ** \param[in]  dlim     Maximum distance (TEST if not defined)
 **
 ** \param[out]  dtab    Array containing the surface of influence
 **                      (Dimension = Number of samples in db_point)
 ** \param[out]  gtab    Array containing the surface of influence of the
 **                      polygon to which it belongs (or TEST)
 **                      (Dimension = Number of samples in db_grid)
 **
 *****************************************************************************/
int surface(Db *db_point,
            DbGrid *db_grid,
            int /*icol*/,
            double dlim,
            double *dtab,
            double *gtab)
{
  bool flagTest;

  if (!db_grid->hasSameDimension(db_point)) return (1);
  int ndim = db_point->getNDim();

  /* Preliminary calculations */

  double d2max = (FFFF(dlim)) ? 1.e30 : dlim * dlim;
  double maille = db_grid->getCellSize();
  for (int iech = 0; iech < db_point->getSampleNumber(); iech++)
    dtab[iech] = 0.;

  /* Loop on the target points */

  VectorDouble vgrid(ndim);
  for (int igrid = 0; igrid < db_grid->getSampleNumber(); igrid++)
  {
    gtab[igrid] = -1;
    if (!db_grid->isActive(igrid)) continue;
    flagTest = false;
    for (int idim = 0; idim < ndim && ! flagTest; idim++)
    {
      vgrid[idim] = db_grid->getCoordinate(igrid, idim);
      if (FFFF(vgrid[idim])) flagTest = true;
    }
    if (flagTest) continue;

    /* Loop on the data points */

    double d2min = d2max;
    for (int iech = 0; iech < db_point->getSampleNumber(); iech++)
    {
      if (!db_point->isActive(iech)) continue;

      /* Calculate the distance between node and data */

      double dist = 0.;
      flagTest = false;
      for (int idim = 0; idim < ndim && ! flagTest; idim++)
      {
        double v2 = db_point->getCoordinate(iech, idim);
        if (FFFF(v2))
          flagTest = true;
        else
        {
          double delta = v2 - vgrid[idim];
          dist += delta * delta;
        }
      }
      if (flagTest) continue;
      if (dist > d2max) continue;

      /* Keep the closest sample */

      if (dist > d2min) continue;
      gtab[igrid] = iech;
      d2min = dist;
    }
  }

  /* Calculate the influence of each datum */

  for (int igrid = 0; igrid < db_grid->getSampleNumber(); igrid++)
  {
    int jech = (int) gtab[igrid];
    if (jech >= 0) dtab[jech]++;
  }
  for (int iech = 0; iech < db_point->getSampleNumber(); iech++)
    dtab[iech] *= maille;

  /* Evaluate each grid node with the size of the influence polygon */
  /* to which it belongs                                            */

  for (int igrid = 0; igrid < db_grid->getSampleNumber(); igrid++)
  {
    int jech = (int) gtab[igrid];
    if (jech >= 0)
      gtab[igrid] = dtab[jech];
    else
      gtab[igrid] = TEST;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Print the Db contents
 **
 ** \param[in]  db    Db descriptor
 ** \param[in]  nrdv  Radius of the variable display
 ** \param[in]  nrds  Radius of the sample display
 ** \param[in]  ivar  Rank of the Target Variable
 ** \param[in]  iech  Rank of the Target Sample
 **
 *****************************************************************************/
static void st_edit_display(Db *db, int nrdv, int nrds, int ivar, int iech)
{
  int item, nvar, nech, ivar_deb, ivar_fin, iech_deb, iech_fin, jvar, jech;
  ELoc locatorType;
  char string[5];

  /* Initializations */

  (void) gslStrcpy(string, "NA");
  nech = db->getSampleNumber();
  nvar = db->getColumnNumber();

  ivar_deb = ivar - nrdv;
  ivar_fin = ivar + nrdv;
  if (ivar_deb < 0)
  {
    ivar_deb = 0;
    ivar_fin = MIN(2 * nrdv, nvar - 1);
  }
  if (ivar_fin >= nvar)
  {
    ivar_fin = nvar - 1;
    ivar_deb = MAX(0, ivar_fin - 2 * nrdv);
  }

  iech_deb = iech - nrds;
  iech_fin = iech + nrds;
  if (iech_deb < 0)
  {
    iech_deb = 0;
    iech_fin = MIN(2 * nrds, nech - 1);
  }
  if (iech_fin >= nech)
  {
    iech_fin = nech - 1;
    iech_deb = MAX(0, iech_fin - 2 * nrds);
  }

  /* Print the Header (Variable name) */

  tab_prints(NULL, " ");
  for (jvar = ivar_deb; jvar <= ivar_fin; jvar++)
  {
    if (db->getLocatorByColIdx(jvar, &locatorType, &item))
    {
      String strloc = getLocatorName(locatorType, item);
      (void) gslStrcpy(string, strloc.c_str());
    }
    else
      (void) gslStrcpy(string, "NA");
    if (jvar == ivar) (void) gslStrcat(string, "*");
    tab_prints(NULL, string);
  }
  message("\n");

  /* Print the Header (Variable rank) */

  tab_prints(NULL, " ");
  for (jvar = ivar_deb; jvar <= ivar_fin; jvar++)
    tab_print_rc(NULL, 2, jvar + 1);
  message("\n");

  /* Loop on the samples */

  for (jech = iech_deb; jech <= iech_fin; jech++)
  {
    tab_print_rc(NULL, 3, jech + 1);
    if (iech == jech)
      message("*");
    else
      message(" ");
    for (jvar = ivar_deb; jvar <= ivar_fin; jvar++)
      tab_printg(NULL, db->getArray(jech, jvar));
    message("\n");
  }
}

/****************************************************************************/
/*!
 **  Look for the next sample with a value within the interval
 **
 ** \return  Rank of the next sample
 **
 ** \param[in]  db      Db descriptor
 ** \param[in]  iech    Rank of the current sample
 ** \param[in]  ivar    Rank of the current variable
 ** \param[in]  orient  Orientation (1: downwards; -1: forwards)
 ** \param[in]  vmin    Minimum value
 ** \param[in]  vmax    Maximum value
 **
 *****************************************************************************/
static int st_edit_find(Db *db,
                        int iech,
                        int ivar,
                        int orient,
                        double vmin,
                        double vmax)
{
  double value;

  /* Dispatch */

  if (orient > 0)
  {
    for (int i = iech + 1; i < db->getSampleNumber(); i++)
    {
      value = db->getArray(i, ivar);
      if (FFFF(value)) continue;
      if (!FFFF(vmin) && value < vmin) continue;
      if (!FFFF(vmax) && value > vmax) continue;
      return (i);
    }
    messerr("--> String not found before the end-of-file");
    return (iech);
  }
  for (int i = iech - 1; i >= 0; i--)
  {
    value = db->getArray(i, ivar);
    if (FFFF(value)) continue;
    if (!FFFF(vmin) && value < vmin) continue;
    if (!FFFF(vmax) && value > vmax) continue;
    return (i);
  }
  messerr("--> String not found before the top-of-file");
  return (iech);
}

/****************************************************************************/
/*!
 **  Ask for the next keyword in the Editor
 **
 ** \return  Return code:
 ** \return  0 : A valid keyword has been found
 ** \return  1 : The 'stop' has been met
 ** \return -1 : The 'quit' has been met
 **
 ** \param[out]  item   Selected item
 ** \param[out]  rank   Value for the Shift
 ** \param[out]  vmin   Value for the lower bound
 ** \param[out]  vmax   Value for the upper bound
 **
 *****************************************************************************/
static int st_edit_ask(int *item, int *rank, double *vmin, double *vmax)
{
  int found, flag_skip, mem_long;
  char string[STRING_LENGTH], *decode;
  static int mem_item = 1;
  static int mem_rank = 1;
  static double mem_vmin = 0.;
  static double mem_vmax = 1.;

  label_loop: _lire_string("Enter Command (or 'stop' or 'quit' or '?')", 0,
  NULL,
                           string);

  /* Look for the string */

  found = -1;
  for (int i = 0; i < N_EDIT; i++)
    if (!strncasecmp(string, EDIT[i].key, strlen(EDIT[i].key))) found = i;

  /* Check for the special keyword */

  if (!strcasecmp(string, "STOP")) return (1);
  if (!strcasecmp(string, "QUIT")) return (-1);

  /* A valid keyword has not been found */

  if (found < 0)
  {
    mestitle(1, "List of the Valid Editor Keywords:");
    for (int i = 0; i < N_EDIT; i++)
      message("%2s : %s\n", EDIT[i].key, EDIT[i].title);
    goto label_loop;
  }

  /* A valid keyword has been encountered: Check for the rest */

  decode = &string[strlen(EDIT[found].key)];

  /* The 'same' command has been encountered */

  flag_skip = 0;
  if (found == 1)
  {
    found = mem_item;
    *rank = mem_rank;
    *vmin = mem_vmin;
    *vmax = mem_vmax;
    flag_skip = 1;
  }

  /* Ask for complementary information */

  if (!flag_skip)
  {
    /* A rank must be specified */

    if (EDIT[found].flag_rank)
    {
      string_strip_blanks(decode, 1);
      mem_long = static_cast<int>(strlen(decode));
      if (mem_long > 0)
      {
        *rank = strtol(decode, &decode, 0);
        if (mem_long == (int) strlen(decode))
        {
          messerr("Cannot convert '%s' into a valid Rank", decode);
          goto label_loop;
        }
      }
      else
        *rank = _lire_int("Value for the Shift", 1, mem_rank, ITEST, ITEST);
    }

    /* Bounds must be specified */

    if (EDIT[found].flag_bounds)
    {
      string_strip_blanks(decode, 1);
      mem_long = static_cast<int>(strlen(decode));
      if (mem_long > 0)
      {
        *vmin = strtod(decode, &decode);
        if (mem_long == (int) strlen(decode))
        {
          messerr("Cannot convert '%s' into a valid Minimum Bound", decode);
          goto label_loop;
        }
      }
      else
        *vmin = _lire_double("Minimum value", 1, mem_vmin, TEST, TEST);

      string_strip_blanks(decode, 1);
      mem_long = static_cast<int>(strlen(decode));
      if (mem_long > 0)
      {
        *vmax = strtod(decode, &decode);
        if (mem_long == (int) strlen(decode))
        {
          messerr("Cannot convert '%s' into a valid Maximum Bound", decode);
          goto label_loop;
        }
        if (*vmax < *vmin)
        {
          messerr("Upper bound (%lf) may not be smaller than Lower bound (%lf)",
                  (*vmax), (*vmin));
        }
      }
      else
        *vmax = _lire_double("Maximum value", 1, mem_vmax, *vmin, TEST);
    }
  }

  /* Return argument */

  *item = found;

  /* Store the answers for default values in next operation */

  mem_item = *item;
  mem_rank = *rank;
  mem_vmin = *vmin;
  mem_vmax = *vmax;

  return (0);
}

/****************************************************************************/
/*!
 **  Edit the Data Base Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db   Db descriptor
 **
 ** \param[out] flag_valid: 1 for 'stop' and 0 for 'quit'
 **
 *****************************************************************************/
int db_edit(Db *db, int *flag_valid)

{
  int nech, nvar, ivar, iech, incr, type, ok, nrds, nrdv, flag_inter;
  double vmin, vmax, value;

  /* Initializations */

  nech = db->getSampleNumber();
  nvar = db->getColumnNumber();
  ivar = iech = 0;
  nrds = nrdv = incr = 1;
  vmin = vmax = TEST;
  if (nech < 1 || nvar < 1) return (1);

  ok = 1;
  while (ok)
  {
    st_edit_display(db, nrdv, nrds, ivar, iech);
    flag_inter = st_edit_ask(&type, &incr, &vmin, &vmax);
    if (flag_inter > 0)
    {
      *flag_valid = 1;
      ok = 0;
      break;
    }
    if (flag_inter < 0)
    {
      *flag_valid = 0;
      ok = 0;
      break;
    }

    /* Dispatch */

    switch (type)
    {
      case 0: /* Set the Parameters */
        nrdv = _lire_int("Display Radius along Variable", 1, nrdv, 0, ITEST);
        nrds = _lire_int("Display Radius along Sample", 1, nrds, 0, ITEST);
        break;

      case 2: /* Relative Sample Rank */
        iech = MAX(0, MIN(iech + incr, nech - 1));
        break;

      case 3: /* Relative Variable Rank */
        ivar = MAX(0, MIN(ivar + incr, nvar - 1));
        break;

      case 4: /* Absolute Sample Rank */
        iech = MAX(0, MIN(incr - 1, nech - 1));
        break;

      case 5: /* Absolute Variable Rank */
        ivar = MAX(0, MIN(incr - 1, nvar - 1));
        break;

      case 6: /* Modify the Value */
        value = _lire_double("New value", 1, db->getArray(iech, ivar), TEST,
        TEST);
        db->setArray(iech, ivar, value);
        break;

      case 7: /* Next sample within an interval */
        iech = st_edit_find(db, iech, ivar, 1, vmin, vmax);
        break;

      case 8: /* Previous sample within an interval */
        iech = st_edit_find(db, iech, ivar, -1, vmin, vmax);
        break;

      case 9: /* Display the current selection */
      default:
        break;
    }
  }
  return (0);
}

/*****************************************************************************/
/*!
 **  Generates the discretized points along the trace
 **
 ** \param[in]  nseg   Number of vertices along the trace
 ** \param[in]  trace  Array defining the trace
 **                    (Dimension: 2 * nseg)
 ** \param[in]  disc   Discretization distance
 **
 ** \param[out] np_arg   Number of discretized points
 ** \param[out] xp_arg   Array of first coordinates
 ** \param[out] yp_arg   Array of second coordinates
 ** \param[out] dd_arg   Array of distances between discretized points
 ** \param[out] del_arg  Array of distances between vertices
 ** \param[out] dist_arg Total distance of the trace
 **
 *****************************************************************************/
void ut_trace_discretize(int nseg,
                         const double *trace,
                         double disc,
                         int *np_arg,
                         double **xp_arg,
                         double **yp_arg,
                         double **dd_arg,
                         double **del_arg,
                         double *dist_arg)
{
  double *xp, *yp, *dd, *del, deltax, deltay, x0, y0, x1, y1, dist;
  int iseg, np, ecr, nloc, ip;

  /* Initializations */

  xp = yp = dd = nullptr;
  (*np_arg) = np = 0;
  (*dist_arg) = x1 = y1 = 0.;
  del = (double*) mem_alloc(sizeof(double) * nseg, 1);
  del[0] = 0.;

  /* Loop on the trace segments */

  for (iseg = ecr = 0; iseg < nseg - 1; iseg++)
  {

    /* Consider a segment trace */

    x0 = TRACE(0, iseg);
    y0 = TRACE(1, iseg);
    x1 = TRACE(0, iseg + 1);
    y1 = TRACE(1, iseg + 1);
    deltax = x1 - x0;
    deltay = y1 - y0;
    dist = sqrt(deltax * deltax + deltay * deltay);
    (*dist_arg) += dist;
    del[iseg + 1] = (*dist_arg);

    /* Discretize the trace segment */

    nloc = (int) floor(dist / disc);
    if (ABS(nloc * disc - dist) < dist / 1000) nloc--;
    np += nloc;
    xp = (double*) mem_realloc((char* ) xp, sizeof(double) * np, 1);
    yp = (double*) mem_realloc((char* ) yp, sizeof(double) * np, 1);

    for (ip = 0; ip < nloc; ip++, ecr++)
    {
      xp[ecr] = x0 + deltax * ip / nloc;
      yp[ecr] = y0 + deltay * ip / nloc;
    }
  }

  /* Adding the last vertex */

  np++;
  xp = (double*) mem_realloc((char* ) xp, sizeof(double) * np, 1);
  yp = (double*) mem_realloc((char* ) yp, sizeof(double) * np, 1);
  xp[ecr] = x1;
  yp[ecr] = y1;
  ecr++;

  /* Elaborate the vector of distances */

  dd = (double*) mem_alloc(sizeof(double) * np, 1);
  dd[0] = 0.;
  for (ip = 0; ip < np - 1; ip++)
  {
    deltax = xp[ip + 1] - xp[ip];
    deltay = yp[ip + 1] - yp[ip];
    dd[ip + 1] = dd[ip] + sqrt(deltax * deltax + deltay * deltay);
  }

  /* Returning arguments */

  (*np_arg) = np;
  (*xp_arg) = xp;
  (*yp_arg) = yp;
  (*dd_arg) = dd;
  (*del_arg) = del;
}

/*****************************************************************************/
/*!
 **  Sample the point Db close to discretized points of the trace
 **
 ** \param[in]  db     Db to be sampled
 ** \param[in]  ptype  Type of locator
 ** \param[in]  np     Number of discretized points
 ** \param[in]  xp     Array of first coordinates
 ** \param[in]  yp     Array of second coordinates
 ** \param[in]  dd     Array of distances
 ** \param[in]  radius Neighborhood radius
 **
 ** \param[out] ns_arg     Number of sampled points
 ** \param[out] xs_arg     Array of first coordinates of sampled points
 ** \param[out] ys_arg     Array of second coordinates of sampled points
 ** \param[out] rks_arg    Array of sample indices (starting from 1)
 ** \param[out] lys_arg    Array of layer indices of sampled points
 ** \param[out] typ_arg    Array of data type
 **                        1 for hard data in Z or TIME
 **                        2 for lower bound
 **                        3 for upper bound
 **
 *****************************************************************************/
void ut_trace_sample(Db *db,
                     const ELoc& ptype,
                     int np,
                     const double *xp,
                     const double *yp,
                     const double *dd,
                     double radius,
                     int *ns_arg,
                     double **xs_arg,
                     double **ys_arg,
                     int **rks_arg,
                     int **lys_arg,
                     int **typ_arg)
{
  int *lys, *typ, *rks, iech, ip, ns, ipmin, nvar;
  double *xs, *ys, cote, layer, bound[2];
  double radcarre, xx, yy, delx, dely, dist, ddmin;

  /* Initializations */

  radcarre = radius * radius;
  xs = ys = nullptr;
  lys = typ = rks = nullptr;
  ns = 0;
  nvar = db->getIntervalNumber();

  /* Loop on the samples */

  for (iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Coordinates of the sample point */

    xx = db->getCoordinate(iech, 0);
    yy = db->getCoordinate(iech, 1);

    /* Loop on the discretized samples */

    ipmin = -1;
    ddmin = 1.e30;
    for (ip = 0; ip < np; ip++)
    {
      delx = xx - xp[ip];
      dely = yy - yp[ip];
      dist = delx * delx + dely * dely;
      if (dist > radcarre || dist > ddmin) continue;
      ddmin = dist;
      ipmin = ip;
    }
    if (ipmin < 0) continue;

    /* Keep sample defined by locator */

    cote = db->getFromLocator(ptype, iech);
    if (!FFFF(cote))
    {
      layer = db->getFromLocator(ELoc::LAYER, iech);
      xs = (double*) mem_realloc((char* ) xs, sizeof(double) * (ns + 1), 1);
      ys = (double*) mem_realloc((char* ) ys, sizeof(double) * (ns + 1), 1);
      lys = (int*) mem_realloc((char* ) lys, sizeof(int) * (ns + 1), 1);
      typ = (int*) mem_realloc((char* ) typ, sizeof(int) * (ns + 1), 1);
      rks = (int*) mem_realloc((char* ) rks, sizeof(int) * (ns + 1), 1);
      xs[ns] = dd[ipmin];
      ys[ns] = cote;
      lys[ns] = (FFFF(layer)) ? 1 : (int) layer + 1;
      typ[ns] = 1;
      rks[ns] = iech + 1;
      ns++;
    }

    /* Keep sample defined by locator UP or LOW */

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      bound[0] = db->getLocVariable(ELoc::L,iech, ivar);
      bound[1] = db->getLocVariable(ELoc::U,iech, ivar);
      for (int ib = 0; ib < 2; ib++)
      {
        if (FFFF(bound[ib])) continue;
        xs = (double*) mem_realloc((char* )xs, sizeof(double) * (ns + 1), 1);
        ys = (double*) mem_realloc((char* )ys, sizeof(double) * (ns + 1), 1);
        lys = (int*) mem_realloc((char* )lys, sizeof(int) * (ns + 1), 1);
        typ = (int*) mem_realloc((char* )typ, sizeof(int) * (ns + 1), 1);
        rks = (int*) mem_realloc((char* )rks, sizeof(int) * (ns + 1), 1);
        xs[ns] = dd[ipmin];
        ys[ns] = bound[ib];
        lys[ns] = ivar + 1;
        typ[ns] = ib + 2;
        rks[ns] = iech + 1;
        ns++;
      }
    }
  }

  /* Returning arguments */

  *ns_arg = ns;
  *xs_arg = xs;
  *ys_arg = ys;
  *lys_arg = lys;
  *typ_arg = typ;
  *rks_arg = rks;
}

/*****************************************************************************/
/*!
 **  Create a set of samples according to a Poisson process
 **
 ** \return  Array of returned values (or NULL)
 **
 ** \param[in]  number      Number of samples to be generated
 ** \param[in]  coormin     Vector of minimum coordinates
 ** \param[in]  coormax     Vector of maximum coordinates
 ** \param[in]  flag_repulsion True if the repulsion process is active
 ** \param[in]  range       Repulsion range
 ** \param[in]  beta        Repulsion beta coefficient
 **
 *****************************************************************************/
static VectorDouble st_point_init_homogeneous(int number,
                                              const VectorDouble &coormin,
                                              const VectorDouble &coormax,
                                              bool flag_repulsion,
                                              double range,
                                              double beta)
{
  VectorDouble tab;

  if (coormin.empty() || coormax.empty())
  {
    messerr("This method requires 'coormin' and 'coormax' defined");
    return tab;
  }
  VectorDouble extend = VH::subtract(coormin, coormax);
  int ndim = (int) coormin.size();
  VectorDouble coor(ndim);
  VectorDouble delta(ndim);

  // Generate the samples

  tab.resize(ndim * number);

  int ecr = 0;
  for (int ip = 0; ip < number; ip++)
  {
    for (int idim = 0; idim < ndim; idim++)
      coor[idim] = coormin[idim] + law_uniform(0., extend[idim]);

    // Check if the point is acceptable

    bool flag_drop = false;
    if (flag_repulsion)
    {

      // Calculate the shortest distance with the previous samples

      double ddmin = 1.e30;
      for (int jp = 0; jp < ip; jp++)
      {
        for (int idim = 0; idim < ndim; idim++)
          delta[idim] = (tab[ndim * jp + idim] - coor[idim]) / range;
        double dd = VH::norm(delta);
        if (dd < ddmin) ddmin = dd;
      }

      /* Check the rejection criterion */

      double proba = exp(-pow(ddmin, beta));
      double alea = law_uniform(0., 1.);
      flag_drop = (alea < proba);
    }

    // Add the new point
    if (flag_drop) continue;
    for (int idim = 0; idim < ndim; idim++)
      tab[ndim * ecr + idim] = coor[idim];
    ecr++;
  }
  tab.resize(ndim * ecr);
  return (tab);
}

/*****************************************************************************/
/*!
 **  Create a set of samples according to a Poisson Regionalized process
 **
 ** \return  Array of returned values (or NULL)
 **
 ** \param[in]  number      Number of samples to be generated
 ** \param[in]  dbgrid      Descriptor of the Db grid parameters
 ** \param[in]  flag_repulsion True if the repulsion process is active
 ** \param[in]  range       Repulsion range
 ** \param[in]  beta        Repulsion beta coefficient
 **
 ** \remarks Thinning can only be defined in 2-D.
 ** \remarks If the thinning is regionalized, its parameters are stored
 ** \remarks  as NOSTAT variables: Range-1, Range-2 and Angle
 **
 *****************************************************************************/
static VectorDouble st_point_init_inhomogeneous(int number,
                                                DbGrid *dbgrid,
                                                bool flag_repulsion,
                                                double range,
                                                double beta)
{
  VectorDouble tab;

  int ndim = dbgrid->getNDim();
  if (dbgrid == nullptr)
  {
    messerr("This function requires a DbGrid data base");
    return tab;
  }
  if (! dbgrid->isGrid())
  {
    messerr("This function requires the Db organized as a grid");
    return tab;
  }
  bool flag_dens = (dbgrid->getLocNumber(ELoc::Z) == 1);
  bool flag_region = (ndim == 2 && dbgrid->getLocatorNumber(ELoc::NOSTAT) == (ndim+1));

  VectorDouble coor(ndim);
  VectorDouble coorbis(ndim);
  VectorDouble delta(ndim);
  VectorDouble radius(ndim);
  VectorDouble radip(ndim);
  double angip = 0.;

  /* Evaluate the density */

  int ngrid = dbgrid->getSampleNumber(true);
  VectorDouble dens;
  dens.resize(ngrid,0.);
  double denstot = 0.;
  if (flag_dens)
  {
    int ig = 0;
    for (int jg = 0, ng = dbgrid->getActiveSampleNumber(); jg < ng; jg++)
    {
      if (!dbgrid->isActiveAndDefined(jg, 0)) continue;
      double densloc = dbgrid->getZVariable(jg, 0);
      if (densloc >= 0) denstot += densloc;
      dens[ig++] = denstot;
    }
  }
  else
  {
    denstot = dbgrid->getSampleNumber(true);
  }

  /* Point generation */

  int ecr = 0;
  int indip = 0;
  int indjp = 0;
  int ntrial = 0;
  while (number - ecr > ntrial / 10)
  {
    // Draw a probability

    double proba = law_uniform(0., denstot);
    ntrial++;

    // Draw a grid cell at random

    if (flag_dens)
    {
      double denscum = 0.;
      indip = -1;
      for (int ig = 0; ig < ngrid && indip < 0; ig++)
      {
        if (!dbgrid->isActive(ig)) continue;
        double densloc = dbgrid->getZVariable(ig, 0);
        if (FFFF(densloc) || densloc < 0) continue;
        denscum += densloc;
        if (denscum > proba) indip = ig;
      }
      if (indip < 0) indip = ngrid - 1;
    }
    else
    {
      indip = (int) proba;
    }

    // Draw the point within the elected cell

    dbgrid->rankToCoordinatesInPlace(indip, coor);
    for (int idim = 0; idim < ndim; idim++)
      coor[idim] += law_uniform(0., dbgrid->getDX(idim));
    if (flag_region)
    {
      for (int idim = 0; idim < ndim; idim++)
        radip[idim] = dbgrid->getFromLocator(ELoc::NOSTAT,indip,idim);
      angip = dbgrid->getFromLocator(ELoc::NOSTAT, indip, ndim);
    }

    // Check if the point is acceptable

    bool flag_drop = false;
    if (flag_repulsion)
    {

      // Calculate the shortest distance with the previous samples

      flag_drop = false;
      for (int jp = 0; jp < ecr && ! flag_drop; jp++)
      {
        double dd = 0.;
        for (int idim = 0; idim < ndim; idim++)
        {
          coorbis[idim] = tab[ndim * jp + idim];
          delta[idim] = (coorbis[idim] - coor[idim]);
        }

        if (! flag_region)
        {
          dd = VH::norm(delta) / range;
        }
        else
        {
          indjp = dbgrid->coordinateToRank(coorbis);
          for (int idim = 0; idim < ndim; idim++)
            radius[idim] = 2. / (radip[idim] + dbgrid->getFromLocator(ELoc::NOSTAT,indjp,idim));
          double angle = (angip + dbgrid->getFromLocator(ELoc::NOSTAT, indjp, ndim)) / 2.;
          Tensor tensor(ndim);
          tensor.setRotationAngle(0,angle);
          tensor.setRadiusVec(radius);
          dd = VH::norm(tensor.applyInverse(delta));
        }

        // Check if the point 'ip' must be dropped
        proba = exp(-pow(dd, beta));
        flag_drop = (law_uniform(0.,1.) < proba);
      }
    }
    if (flag_drop) continue;

    // Add the new point
    for (int idim = 0; idim < ndim; idim++)
      tab.push_back(coor[idim]);
    ecr++;
  }

  return tab;
}

/****************************************************************************/
/*!
 **  Create indicator residual variables
 **
 ** \return  Error returned code
 **
 ** \param[in]  db      Db structure
 ** \param[in]  ivar    Index of the target variable
 ** \param[in]  zcut    Array containing the cutoffs
 **
 ** \remarks The array 'zcut' must be provided in increasing order
 **
 *****************************************************************************/
int db_resind(Db *db, int ivar, const VectorDouble& zcut)
{
  int nech = db->getSampleNumber();
  int ncut = (int) zcut.size();
  if (! VH::isSorted(zcut, true))
  {
    messerr("The cutoffs must be provided in increasing order");
    return 1;
  }

  /* Calculate the tonnages */

  int ntot = 0;
  VectorDouble tonnage(ncut, 0);
  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    double value = db->getArray(iech, ivar);
    if (FFFF(value)) continue;
    ntot++;

    for (int icut = 0; icut < ncut; icut++)
      if (value >= zcut[icut]) tonnage[icut]++;
  }
  for (int icut = 0; icut < ncut; icut++)
    tonnage[icut] /= (double) ntot;

  /* Create the variables */

  int iptr = db->addColumnsByConstant(ncut, TEST);
  if (iptr < 0) return 1;

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    double value = db->getArray(iech, ivar);
    if (FFFF(value)) continue;

    /* Loop on the cutoffs */

    for (int icut = 0; icut < ncut; icut++)
    {
      double zval = zcut[icut];
      int ind_cut0 = (value > zval);
      zval = (icut > 0) ? zcut[icut - 1] : 0.;
      int ind_cut1 = (value > zval);
      double ton_cut0 = tonnage[icut];
      double ton_cut1 = (icut > 0) ? tonnage[icut - 1] : 1.;
      int ir = ind_cut0 / ton_cut0 - ind_cut1 / ton_cut1;
      db->setArray(iech, iptr + icut, ir);
    }
  }

  return (0);
}

/****************************************************************************/
/*!
 **  Normalize the gradient components
 **
 ** \param[in]  dbgrid  Db structure (grid organized)
 **
 *****************************************************************************/
static void st_gradient_normalize(Db *dbgrid)

{
  double norme, grad;
  int ndim;

  /* Initializations */

  ndim = dbgrid->getNDim();

  /* Loop on the samples */

  for (int iech = 0; iech < dbgrid->getSampleNumber(); iech++)
  {

    norme = 0.;
    for (int idim = 0; idim < ndim; idim++)
    {
      grad = dbgrid->getLocVariable(ELoc::G,iech, idim);
      norme += grad * grad;
    }

    if (norme <= 0) continue;
    norme = sqrt(norme);

    for (int idim = 0; idim < ndim; idim++)
    {
      grad = dbgrid->getLocVariable(ELoc::G,iech, idim);
      dbgrid->setLocVariable(ELoc::G,iech, idim, grad / norme);
    }
  }
}

/****************************************************************************/
/*!
 **  Calculate the gradient over a grid
 **
 ** \return  Rank of the newly created variables (or -1)
 **
 ** \param[in]  dbgrid  Db structure (grid organized)
 **
 *****************************************************************************/
int db_gradient_components(DbGrid *dbgrid)

{
  int iptrz, iptr, nx, ny, nz, nmax, error, ndim, j1, j2, number;
  double dinc, v1, v2, delta;
  VectorInt indg;

  /* Preliminary check */

  error = number = 1;
  iptr = -1;
  ndim = dbgrid->getNDim();
  if (! dbgrid->isGrid())
  {
    messerr("The Db should be organized as a Grid");
    goto label_end;
  }
  if (!dbgrid->isVariableNumberComparedTo(1)) goto label_end;
  if (ndim > 3)
  {
    messerr("This function is limited to Space Dimension <= 3");
    goto label_end;
  }

  /* Initializations */

  nx = dbgrid->getNX(0);
  ny = dbgrid->getNX(1);
  nz = dbgrid->getNX(2);
  indg.resize(ndim, 0);

  /* Create the new variable */

  iptrz = dbgrid->getColIdxByLocator(ELoc::Z, 0);
  if (iptrz < 0) goto label_end;
  iptr = dbgrid->addColumnsByConstant(ndim, TEST, String(), ELoc::G);

  /* Calculate the Gradient components */

  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      for (int iz = 0; iz < nz; iz++)
      {
        for (int idim = 0; idim < ndim; idim++)
        {
          nmax = dbgrid->getNX(idim);
          dinc = dbgrid->getDX(idim);

          v1 = v2 = 0.;
          if (idim == 0)
          {
            j1 = (ix + 1 > nmax - 1) ? ix : ix + 1;
            v1 = get_grid_value(dbgrid, iptrz, indg, j1, iy, iz);
            if (FFFF(v1)) continue;
            j2 = (ix - 1 < 0) ? ix : ix - 1;
            v2 = get_grid_value(dbgrid, iptrz, indg, j2, iy, iz);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          if (idim == 1)
          {
            j1 = (iy + 1 > nmax - 1) ? iy : iy + 1;
            v1 = get_grid_value(dbgrid, iptrz, indg, ix, j1, iz);
            if (FFFF(v1)) continue;
            j2 = (iy - 1 < 0) ? iy : iy - 1;
            v2 = get_grid_value(dbgrid, iptrz, indg, ix, j2, iz);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          if (idim == 2)
          {
            j1 = (iz + 1 > nmax - 1) ? iz : iz + 1;
            v1 = get_grid_value(dbgrid, iptrz, indg, ix, iy, j1);
            if (FFFF(v1)) continue;
            j2 = (iz - 1 < 0) ? iz : iz - 1;
            v2 = get_grid_value(dbgrid, iptrz, indg, ix, iy, j2);
            if (FFFF(v2)) continue;
            number = j1 - j2;
          }

          delta = (v1 - v2) / (number * dinc);
          set_grid_value(dbgrid, iptr + idim, indg, ix, iy, iz, delta);
        }
      }

  /* Set the error return code */

  error = 0;

  label_end: if (error)
  {
    (void) db_attribute_del_mult(dbgrid, iptr, ndim);
    iptr = -1;
  }
  return (iptr);
}

/****************************************************************************/
/*!
 **  Check if one (at least) of the gradient components is undefined
 **
 ** \return  1 If one component (at least) is undefined
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  iptr_grad Rank of the first gradient component
 ** \param[in]  iech      Sample rank
 **
 *****************************************************************************/
static int st_is_undefined(Db *dbgrid, int iptr_grad, int iech)
{
  int ndim;

  ndim = dbgrid->getNDim();
  for (int idim = 0; idim < ndim; idim++)
  {
    if (FFFF(dbgrid->getArray(iech, iptr_grad + idim))) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Check if the gradient is null
 **
 ** \return  1 If gradient is null
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  iptr_grad Rank of the first gradient component
 ** \param[in]  iech      Sample rank
 **
 *****************************************************************************/
static int st_is_zero(Db *dbgrid, int iptr_grad, int iech)
{
  double grad, delta;
  int ndim;

  grad = 0.;
  ndim = dbgrid->getNDim();
  for (int idim = 0; idim < ndim; idim++)
  {
    delta = dbgrid->getArray(iech, iptr_grad + idim);
    grad += delta * delta;
  }
  return (grad < 1.e-5);
}

/****************************************************************************/
/*!
 **  Get the next gradient-based data
 **
 ** \return  1 Error return code
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  iptr_grad Rank of the first gradient component
 ** \param[in]  coor      Array of coordinates
 **
 ** \param[out] knd       Absolute index
 ** \param[out] surf      Local value of the surface
 **
 *****************************************************************************/
static int st_get_next(DbGrid *dbgrid,
                       int iptr_grad,
                       VectorDouble &coor,
                       int *knd,
                       double *surf)
{
  int knd_loc;
  double surf_loc;

  knd_loc = dbgrid->coordinateToRank(coor);
  if (knd_loc < 0) return 1;
  if (!dbgrid->isActive(knd_loc)) return 1;
  surf_loc = dbgrid->getZVariable(knd_loc, 0);
  if (FFFF(surf_loc) || st_is_undefined(dbgrid, iptr_grad, knd_loc)) return (1);
  if (st_is_zero(dbgrid, iptr_grad, knd_loc)) return (1);
  *knd = knd_loc;
  *surf = surf_loc;
  return 0;
}

/****************************************************************************/
/*!
 **  Calculate the streamlines
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid    Db structure (grid organized)
 ** \param[in]  dbpoint   Db structure for control points
 ** \param[in]  niter     Maximum number of iterations
 ** \param[in]  step      Progress step
 ** \param[in]  flag_norm 1 if the gradients must be normalized
 ** \param[in]  use_grad  1 if the existing gradients must be used
 **                       0 the gradients must be calculated here
 ** \param[in]  save_grad 1 save the gradients generated in this function
 **                       0 delete gradients when calculated here
 **
 ** \param[out] nbline_loc Number of streamline steps
 ** \param[out] npline_loc Number of values per line vertex
 ** \param[out] line_loc   Array of streamline steps (Dimension: 5 * nbline)
 **
 ** \remarks The returned array 'line_loc' must be freed by the calling function
 ** \remarks Use get_keypone("Streamline_Skip",1) to define the skipping ratio
 **
 *****************************************************************************/
int db_streamline(DbGrid *dbgrid,
                  Db *dbpoint,
                  int niter,
                  double step,
                  int flag_norm,
                  int use_grad,
                  int save_grad,
                  int *nbline_loc,
                  int *npline_loc,
                  double **line_loc)
{
  int error, npline, idim, ecr;
  int iptr_time, iptr_accu, iptr_grad, nbline, knd, nquant, nbyech, ndim;
  double *coor0, *line, surf, date;
  static int quant = 1000;
  VectorDouble coor;

  /* Initializations */

  error = 1;
  nbline = nquant = 0;
  iptr_grad = -1;
  coor0 = line = nullptr;
  if (dbpoint == nullptr) dbpoint = dbgrid;
  nbyech = (int) get_keypone("Streamline_Skip", 1.);

  /* Preliminary checks */

  ndim = dbgrid->getNDim();
  if (ndim < 2 || ndim > 3)
  {
    messerr("This function is limited to 2-D or 3-D case");
    goto label_end;
  }
  npline = ndim + 3;

  /* Core allocation on the Grid Db */

  coor.resize(ndim);
  coor0 = db_sample_alloc(dbgrid, ELoc::X);
  if (coor0 == nullptr) goto label_end;
  iptr_time = dbgrid->addColumnsByConstant(1, TEST);
  if (iptr_time < 0) goto label_end;
  iptr_accu = dbgrid->addColumnsByConstant(1, 0.);
  if (iptr_accu < 0) goto label_end;

  /* Calculate the gradient */

  if (use_grad)
  {
    if (dbgrid->getLocNumber(ELoc::G) != ndim)
    {
      messerr("When using the option 'use.grad'");
      messerr("the number of gradients should be %d (%d)", ndim,
              dbgrid->getLocNumber(ELoc::G));
      goto label_end;
    }
    iptr_grad = dbgrid->getColIdxByLocator(ELoc::G, 0);
  }
  else
  {
    iptr_grad = db_gradient_components(dbgrid);
  }
  if (iptr_grad < 0) goto label_end;

  /* Normalize the gradient (optional) */

  if (flag_norm) st_gradient_normalize(dbgrid);

  /* Loop on the drop points */

  for (int iech = 0; iech < dbpoint->getSampleNumber(); iech++)
  {
    if (!dbpoint->isActive(iech)) continue;
    if (iech % nbyech != 0) continue;
    db_sample_load(dbpoint, ELoc::X, iech, coor.data());
    if (st_get_next(dbgrid, iptr_grad, coor, &knd, &surf)) break;

    /* Store the new point in the Streamline */

    if (nbline >= nquant * quant)
    {
      nquant++;
      line = (double*) mem_realloc((char* ) line,
                                   sizeof(double) * npline * nquant * quant, 1);
    }
    for (idim = ecr = 0; idim < ndim; idim++)
      LINE(nbline,ecr++) = coor[idim];
    LINE(nbline,ecr++) = surf;
    LINE(nbline,ecr++) = knd + 1.;
    LINE(nbline,ecr++) = 0.;
    nbline++;

    for (int i = 0; i < niter; i++)
    {
      for (idim = 0; idim < ndim; idim++)
        coor[idim] -= step * dbgrid->getArray(knd, iptr_grad + idim);
      if (st_get_next(dbgrid, iptr_grad, coor, &knd, &surf)) break;

      /* Store the new point in the Streamline */

      if (nbline >= nquant * quant)
      {
        nquant++;
        line = (double*) mem_realloc((char* ) line,
                                     sizeof(double) * npline * nquant * quant,
                                     1);
      }
      for (idim = ecr = 0; idim < ndim; idim++)
        LINE(nbline,ecr++) = coor[idim];
      LINE(nbline,ecr++) = surf;
      LINE(nbline,ecr++) = knd + 1.;
      LINE(nbline,ecr++) = i + 1.;
      nbline++;

      /* Update variables in the grid Db */

      date = MIN(dbgrid->getArray(knd, iptr_time), i + 1.);
      dbgrid->setArray(knd, iptr_time, date);
      dbgrid->updArray(knd, iptr_accu, EOperator::ADD, 1.);
    }

    /* Add the endpoint */

    if (nbline >= nquant * quant)
    {
      nquant++;
      line = (double*) mem_realloc((char* ) line,
                                   sizeof(double) * npline * nquant * quant, 1);
    }
    for (idim = ecr = 0; idim < ndim; idim++)
      LINE(nbline,ecr++) = TEST;
    LINE(nbline,ecr++) = TEST;
    LINE(nbline,ecr++) = TEST;
    LINE(nbline,ecr++) = TEST;
    nbline++;
  }

  /* Final reallocation */

  line = (double*) mem_realloc((char* ) line, sizeof(double) * npline * nbline, 1);

  /* Set the error return code */

  *nbline_loc = nbline;
  *npline_loc = npline;
  *line_loc = line;
  error = 0;

  label_end:
  db_sample_free(coor0);
  if (!use_grad && !save_grad && iptr_grad >= 0)
    db_attribute_del_mult(dbgrid, iptr_grad, ndim);
  return (error);
}

/*****************************************************************************/
/*!
 **  Calculate and store new variables in the Db which contain
 **  the non-stationary Model component
 **
 ** \return  Distance value
 **
 ** \param[in]  db          Db structure
 ** \param[in]  model       Model structure
 ** \param[in]  icov        Rank of the Basic structure
 ** \param[in]  namconv     Naming convention
 **
 ** \remarks This procedure automatically creates several fields:
 ** \remarks ndim fields for storing the ranges
 ** \remarks ndim fields for storing the angles
 ** \remarks 1 field for storing the sill
 **
 *****************************************************************************/
int db_model_nostat(Db *db,
                    Model *model,
                    int icov,
                    const NamingConvention &namconv)
{
  if (icov < 0 || icov >= model->getCovaNumber()) return 1;
  if (!model->isNoStat()) return 0;
  ANoStat *nostat = model->getNoStatModify();

  // The Non-stationary must be defined in the tabulated way
  if (nostat->manageInfo(1, db, nullptr)) return 1;

  /* Create the new variables */

  int ndim = model->getDimensionNumber();
  CovInternal covint(1, -1, 1, -1, ndim, db, db);
  int iptr = db->addColumnsByConstant(2 * ndim + 1, 0.);
  if (iptr < 0) return 1;

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Load the non_stationary parameters */

    covint.setIech1(iech);
    covint.setIech2(iech);
    model->nostatUpdate(&covint);
    CovAniso *cova = model->getCova(icov);

    /* Store the variables */

    int jptr = iptr;
    for (int idim = 0; idim < ndim; idim++)
    {
      db->setArray(iech, jptr, cova->getRange(idim));
      jptr++;
    }
    for (int idim = 0; idim < ndim; idim++)
    {
      db->setArray(iech, jptr, cova->getAnisoAngles(idim));
      jptr++;
    }
    db->setArray(iech, jptr++, cova->getSill(0, 0));
  }

  // Naming convention

  int jptr = iptr;
  for (int idim = 0; idim < ndim; idim++)
    namconv.setNamesAndLocators(
        nullptr, VectorString(), ELoc::UNKNOWN, -1, db, jptr++,
        concatenateStrings("-", "Range", toString(idim + 1)));
  for (int idim = 0; idim < ndim; idim++)
    namconv.setNamesAndLocators(
        nullptr, VectorString(), ELoc::UNKNOWN, -1, db, jptr++,
        concatenateStrings("-", "Angle", toString(idim + 1)));
  namconv.setNamesAndLocators(nullptr, VectorString(), ELoc::UNKNOWN, -1, db, jptr++, "Sill");
  namconv.setLocators(db, iptr, 1, 2 * ndim + 1);

  (void) nostat->manageInfo(-1, db, nullptr);
  return 0;
}

/*****************************************************************************/
/*!
 **  Smooth out the VPC
 **
 ** \return  Distance value
 **
 ** \param[in]  db          3-D Db structure containing the VPCs
 ** \param[in]  width       Width of the Filter
 ** \param[in]  range       Range of the Gaussian Weighting Function
 **
 ** \remarks Work is performed IN PLACE
 **
 *****************************************************************************/
int db_smooth_vpc(DbGrid *db, int width, double range)
{
  int iz, nz, nprop, ecr, nkern, jz, error;
  double *prop1, *prop2, *kernel, total, propval, dz, quant, quant0;

  /* Initializations */

  error = 1;
  nprop = db->getLocNumber(ELoc::P);
  nz = db->getNX(2);
  dz = db->getDX(2);
  prop1 = prop2 = kernel = nullptr;

  /* Core allocation */

  quant0 = law_invcdf_gaussian(0.975);
  if (FFFF(range))
    range = dz * width / quant0;
  else if (IFFFF(width))
    width = (int) (range * quant0 / dz);
  else
  {
    messerr("You must define either 'width' or 'range'");
    goto label_end;
  }
  nkern = 2 * width + 1;
  prop1 = (double*) mem_alloc(sizeof(double) * nz * nprop, 1);
  prop2 = (double*) mem_alloc(sizeof(double) * nz * nprop, 1);
  kernel = (double*) mem_alloc(sizeof(double) * nkern, 1);

  /* Establish the Kernel */

  for (int i = -width; i <= width; i++)
  {
    quant = (i * dz) / range;
    kernel[i + width] = law_df_gaussian(quant) / range;
  }

  /* Preliminary checks */

  if (! db->isGrid() || db->getNDim() != 3) goto label_end;

  /* Loop on the 2-D grid cells */

  for (int ix = 0; ix < db->getNX(0); ix++)
    for (int iy = 0; iy < db->getNX(1); iy++)
    {

      /* Load the proportions */

      if (db_prop_read(db, ix, iy, prop1)) goto label_end;

      /* Loop on the proportions */

      for (int iprop = 0; iprop < nprop; iprop++)
      {

        /* Loop on the levels of the VPC */

        for (iz = ecr = 0; iz < nz; iz++)
        {

          /* Loop on the kernel items */

          total = 0.;
          for (int i = -width; i <= width; i++)
          {
            jz = Grid::generateMirrorIndex(nz, iz+i);
            propval = PROP1(jz, iprop);
            total += kernel[i + width] * propval;
          }
          PROP2(iz,iprop) = total;
        }
      }
      if (db_prop_write(db, ix, iy, prop2)) goto label_end;
    }

  /* Set the error return code */

  error = 0;

  label_end:
  mem_free((char* ) prop1);
  mem_free((char* ) prop2);
  mem_free((char* ) kernel);
  return (error);
}

/*****************************************************************************/
/*!
 **  Regularize variables along vertical wells
 **
 ** \return  Pointer to the newly created Db
 **
 ** \param[in]  db          Initial Db
 ** \param[in]  dbgrid      Reference Grid
 ** \param[in]  flag_center When TRUE, the sample is centered in the layer
 **                         to which it belongs
 **
 ** \remarks This function requires the input well ('db') and the grid to be
 ** \remarks defined in space >= 3D
 ** \remarks It requires a CODE variable to be defined in the input 'db'
 ** \remarks This function regularizes all the variables marked with a Z-locator
 ** \remarks This function takes a sample into account only if isotopic
 **
 *****************************************************************************/
Db* db_regularize(Db *db, DbGrid *dbgrid, int flag_center)
{
  Db* dbnew = nullptr;
  if (db == nullptr) return dbnew;

  // Preliminary checks */

  if (! dbgrid->isGrid())
  {
    messerr("This function requires 'dbgrid' to correspond to a Grid");
    return dbnew;
  }

  if (db->getNDim() < 3)
  {
    messerr("This function requires the 'db' to be defined in 3D or more");
    return dbnew;
  }

  if (dbgrid->getNDim() < 3)
  {
    messerr("This function requires the 'dbgrid' to be defined in 3D or more");
    return dbnew;
  }

  if (!db->hasLocVariable(ELoc::C))
  {
    messerr("This function requires the definition of a CODE variable in 'db'");
    return dbnew;
  }

  if (! db->isVariableNumberComparedTo(1,1))
  {
    messerr("You should define some Z-variables in input 'db'");
    return dbnew;
  }

  // Core allocation 

  int iz = 0;
  int nz   = dbgrid->getNX(2);
  int nvar = db->getLocNumber(ELoc::Z);
  int ndim = db->getNDim();
  int size = ndim + nvar + 1;

  VectorDouble codes = db->getCodeList();
  int ncode = (int) codes.size();
  VectorDouble coor(ndim, 0);
  VectorDouble wcnt(ncode * nz, 0);
  VectorDouble wcor(ncode * nz * ndim, 0);
  VectorDouble wtab(ncode * nz * nvar, 0);

  // Loop on the different samples

  int ntot = db->getSampleNumber();

  //message("Before regularization: ncode = %d, nz = %d, ntot = %d\n", (int)ncode, (int)nz, (int)ntot);

  for (int iech = 0; iech < ntot; iech++)
  {
    if (!db->isActive(iech)) continue;
    mes_process("Regularize Wells", ntot, iech);
    int code = db->getLocVariable(ELoc::C,iech,0);

    // Identify the rank of the code

    int icode = -1;
    for (int i = 0; i < ncode && icode < 0; i++)
      if (areEqual(code,codes[i])) icode = i;
    if (icode < 0) continue;

    // Load the coordinates

    for (int idim = 0; idim < ndim; idim++)
      coor[idim] = db->getCoordinate(iech, idim);

    int err = point_to_bench(dbgrid, coor.data(), 0, &iz);
    if (err < 0) continue;
    if (iz < 0 || iz >= nz) continue;

    // Check if all variables are defined

    int not_defined = 0;
    for (int ivar = 0; ivar < nvar && not_defined == 0; ivar++)
      if (FFFF(db->getZVariable(iech, ivar))) not_defined = 1;
    if (not_defined) continue;

    // Cumulate this sample

    WCNT(iz,icode) += 1.;
    for (int idim = 0; idim < ndim; idim++)
      WCOR(iz,icode,idim) += db->getCoordinate(iech, idim);
    for (int ivar = 0; ivar < nvar; ivar++)
      WTAB(iz,icode,ivar) += db->getZVariable(iech, ivar);
  }

  // Normalization

  int nech = 0;
  for (int icode = 0; icode < ncode; icode++)
    for (iz = 0; iz < nz; iz++)
    {
      double ratio = WCNT(iz, icode);
      if (ratio <= 0) continue;
      for (int idim = 0; idim < ndim; idim++)
        WCOR(iz,icode,idim) /= ratio;
      if (flag_center)
      WCOR(iz,icode,2) = dbgrid->getX0(2) + (0.5 + iz) * dbgrid->getDX(2);
      for (int ivar = 0; ivar < nvar; ivar++)
        WTAB(iz,icode,ivar) /= ratio;
      nech++;
    }

  // Load in storing array

  VectorDouble wecr(size * nech);

  int ecr = 0;
  for (int icode = 0; icode < ncode; icode++)
    for (iz = 0; iz < nz; iz++)
    {
      double ratio = WCNT(iz, icode);
      if (ratio <= 0) continue;
      for (int idim = 0; idim < ndim; idim++)
        wecr[ecr++] = WCOR(iz, icode, idim);
      wecr[ecr++] = codes[icode];
      for (int ivar = 0; ivar < nvar; ivar++)
        wecr[ecr++] = WTAB(iz, icode, ivar);
    }

  // Create the new db

  dbnew = Db::createFromSamples(nech, ELoadBy::SAMPLE, wecr, VectorString(), VectorString(), 0);
  if (dbnew == nullptr) goto label_end;

  ecr = 0;
  dbnew->setLocatorsByUID(ndim, ecr, ELoc::X);
  ecr += ndim;
  dbnew->setLocatorByUID(ecr, ELoc::C);
  ecr += 1;
  dbnew->setLocatorsByUID(nvar, ecr, ELoc::Z);
  ecr += nvar;
  DECLARE_UNUSED(ecr);

  label_end:
  return dbnew;
}

/*****************************************************************************/
/*!
 **  Sampling a fine grid in a coarser set of cells
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbgrid      reference Grid
 ** \param[in]  nvar        Number of variables
 ** \param[in]  vars        Array of variable ranks
 ** \param[in]  npacks      Vector of packing factors
 ** \param[in]  npcell      Number of samples per cell
 ** \param[in]  nmini       Minimum number of nodes before drawing
 **
 ** \param[out] nech_ret    Number of selected samples
 ** \param[out] coor_ret    Array of coordinates
 ** \param[out] data_ret    Array of variables
 **
 ** \remarks The returned arrays 'coor' and 'data' must be freed by
 ** \remarks the calling function
 **
 *****************************************************************************/
int db_grid2point_sampling(DbGrid *dbgrid,
                           int nvar,
                           int *vars,
                           const int *npacks,
                           int npcell,
                           int nmini,
                           int *nech_ret,
                           double **coor_ret,
                           double **data_ret)
{
  int ndim, ntotal, nech, nret, nfine, iech, ecrc, ecrd, error;
  int *retain;
  double *coor, *data;
  VectorInt ranks;
  VectorDouble rndval;

  // Initializations

  *nech_ret = 0;

  error = 1;
  coor = data = nullptr;
  retain = nullptr;
  ndim = dbgrid->getNDim();
  nfine = dbgrid->getSampleNumber();
  nmini       = MAX(nmini, npcell);
  VectorInt indg(ndim,0);
  if (ndim > 3)
  {
    messerr("This function is limited to 3D or less");
    goto label_end;
  }

  // Core allocation 

  ntotal = 1;
  for (int idim = 0; idim < ndim; idim++)
    ntotal *= npacks[idim];
  rndval.resize(ntotal);
  ranks.resize(ntotal);
  retain = (int*) mem_alloc(sizeof(int) * nfine, 0);
  if (retain == nullptr) goto label_end;

  // Dispatch

  nret = 0;
  if (ndim == 1)
  {
    for (int ixcell = 0; ixcell < dbgrid->getNX(0); ixcell += npacks[0])
    {

      // Collect eligible samples

      nech = 0;
      for (int ix = 0; ix < npacks[0]; ix++)
      {
        indg[0] = ixcell + ix;
        if (indg[0] >= dbgrid->getNX(0)) break;
        iech = dbgrid->indiceToRank(indg);
        if (dbgrid->isActive(iech)) ranks[nech++] = iech;
      }
      if (nech < nmini) continue;

      // Draw sample(s) at random 

      for (int i = 0; i < nech; i++)
        rndval[i] = law_uniform(0., 1.);
      VH::arrangeInPlace(0, ranks, rndval, true, nech);
      for (int i = 0; i < npcell; i++)
        retain[nret++] = ranks[i];
    }
  }
  else if (ndim == 2)
  {
    for (int ixcell = 0; ixcell < dbgrid->getNX(0); ixcell += npacks[0])
      for (int iycell = 0; iycell < dbgrid->getNX(1); iycell += npacks[1])
      {

        // Collect eligible samples

        nech = 0;
        for (int ix = 0; ix < npacks[0]; ix++)
          for (int iy = 0; iy < npacks[1]; iy++)
          {
            indg[0] = ixcell + ix;
            if (indg[0] >= dbgrid->getNX(0)) break;
            indg[1] = iycell + iy;
            if (indg[1] >= dbgrid->getNX(1)) break;
            iech = dbgrid->indiceToRank(indg);
            if (dbgrid->isActive(iech)) ranks[nech++] = iech;
          }
        if (nech < nmini) continue;

        // Draw sample(s) at random 

        for (int i = 0; i < nech; i++)
          rndval[i] = law_uniform(0., 1.);
        VH::arrangeInPlace(0, ranks, rndval, true, nech);
        for (int i = 0; i < npcell; i++)
          retain[nret++] = ranks[i];
      }
  }
  else
  {
    for (int ixcell = 0; ixcell < dbgrid->getNX(0); ixcell += npacks[0])
      for (int iycell = 0; iycell < dbgrid->getNX(1); iycell += npacks[1])
        for (int izcell = 0; izcell < dbgrid->getNX(2); izcell += npacks[2])
        {

          // Collect eligible samples

          nech = 0;
          for (int ix = 0; ix < npacks[0]; ix++)
            for (int iy = 0; iy < npacks[1]; iy++)
              for (int iz = 0; iz < npacks[2]; iz++)
              {
                indg[0] = ixcell + ix;
                if (indg[0] >= dbgrid->getNX(0)) break;
                indg[1] = iycell + iy;
                if (indg[1] >= dbgrid->getNX(1)) break;
                indg[2] = izcell + iz;
                if (indg[2] >= dbgrid->getNX(2)) break;
                iech = dbgrid->indiceToRank(indg);
                if (dbgrid->isActive(iech)) ranks[nech++] = iech;
              }
          if (nech < nmini) continue;

          // Draw sample(s) at random 

          for (int i = 0; i < nech; i++)
            rndval[i] = law_uniform(0., 1.);
          VH::arrangeInPlace(0, ranks, rndval, true, nech);
          for (int i = 0; i < npcell; i++)
            retain[nret++] = ranks[i];
        }
  }

  // Allocate the array for coordinates and data

  coor = (double*) mem_alloc(sizeof(double) * ndim * nret, 0);
  if (coor == nullptr) goto label_end;
  data = (double*) mem_alloc(sizeof(double) * nvar * nret, 0);
  if (data == nullptr) goto label_end;

  // Load the returned arrays

  ecrc = ecrd = 0;
  for (int i = 0; i < nret; i++)
  {
    iech = retain[i];
    for (int idim = 0; idim < ndim; idim++)
      coor[ecrc++] = dbgrid->getCoordinate(iech, idim);
    for (int ivar = 0; ivar < nvar; ivar++)
      data[ecrd++] = dbgrid->getArray(iech, vars[ivar]);
  }

  // Set the error return code

  *nech_ret = nret;
  *coor_ret = coor;
  *data_ret = data;
  error = 0;

  // Core deallocation

  label_end:
  retain = (int*) mem_free((char* ) retain);
  return (error);
}

/*****************************************************************************/
/*!
 **  Create a new Data Base with points generated at random
 **
 ** \return  Pointer for the new Db structure
 **
 ** \param[in]  nech        Expected number of samples
 ** \param[in]  coormin     Vector of lower coordinates
 ** \param[in]  coormax     Vector of upper coordinates
 ** \param[in]  dbgrid      Descriptor of the Db grid parameters
 ** \param[in]  flag_exact  True if the number of samples is dran from Poisson
 ** \param[in]  flag_repulsion Use repulsion (need: 'range' and 'beta')
 ** \param[in]  range       Repulsion range
 ** \param[in]  beta        Bending coefficient
 ** \param[in]  extend      Extension of the bounding box (when positive)
 ** \param[in]  seed        Seed for the random number generator
 ** \param[in]  flagAddSampleRank true if the Rank must be generated in the output Db
 **
 ** \remarks Arguments 'extend' is only valid when 'dbgrid' is not defined
 **
 *****************************************************************************/
Db* db_point_init(int nech,
                  const VectorDouble &coormin,
                  const VectorDouble &coormax,
                  DbGrid *dbgrid,
                  bool flag_exact,
                  bool flag_repulsion,
                  double range,
                  double beta,
                  double extend,
                  int seed,
                  bool flagAddSampleRank)
{
  VectorDouble tab;
  Db* db = nullptr;
  int ndim = 0;
  if (dbgrid == nullptr)
    ndim = (int) coormin.size();
  else
    ndim = dbgrid->getNDim();
  if (ndim <= 0) return db;

  // Initiate the pseudo-random number generator

  law_set_random_seed(seed);

  // Process the bounding box extension (optional)

  VectorDouble locmin = coormin;
  VectorDouble locmax = coormax;
  if (extend > 0)
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      locmin[idim] -= extend;
      locmax[idim] += extend;
    }
  }

  // Draw the number of data to be generated in the Poisson process

  int number = nech;
  if (! flag_exact) law_poisson(nech);

  // Dispatch

  if (number > 0)
  {
    if (dbgrid == nullptr)
    {
      tab = st_point_init_homogeneous(number, locmin, locmax,
                                      flag_repulsion, range, beta);
    }
    else
    {
      tab = st_point_init_inhomogeneous(number, dbgrid,
                                        flag_repulsion, range, beta);
    }
  }

  /* Allocate the main structure */

  number = (int) tab.size() / ndim;
  db = Db::createFromSamples(number, ELoadBy::SAMPLE, tab, VectorString(),
                                 VectorString(), flagAddSampleRank);

  /* Set the locators */

  VectorString names = generateMultipleNames("x", ndim);
  for (int idim = 0; idim < ndim; idim++)
  {
    int jdim = (flagAddSampleRank) ? idim + 1 : idim;
    db->setNameByUID(jdim, names[idim]);
    db->setLocatorByUID(jdim, ELoc::X, idim);
  }
  return (db);
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  model       Model structure
 ** \param[in]  niter       Number of iterations
 ** \param[in]  verbose     Verbose flag
 ** \param[in]  namconv     Naming convention
 **
 ** \remarks The procedure uses the FIRST covariance of the Model
 ** \remarks to describe the spatial structure
 **
 *****************************************************************************/
int db_proportion_estimate(Db *dbin,
                           DbGrid *dbout,
                           Model *model,
                           int niter,
                           bool verbose,
                           const NamingConvention &namconv)
{
  VectorVectorInt splits;

  // Preliminary checks

  if (dbin == nullptr)
  {
    messerr("This method requires a 'dbin' argument");
    return 1;
  }
  if (dbout == nullptr)
  {
    messerr("This method requires a 'dbout' argument");
    return 1;
  }
  if (model == nullptr)
  {
    messerr("This method requires a 'model' argument");
    return 1;
  }
  if (dbin->getLocNumber(ELoc::Z) != 1)
  {
    messerr("The argument 'dbin' should have a single variable");
    return 1;
  }

  // Define the environment

  MeshETurbo mesh = MeshETurbo(dbout);
  ShiftOpCs S = ShiftOpCs(&mesh, model, dbout, 0, 0);
  PrecisionOp Qprop = PrecisionOp(&S, model->getCova(0));
  ProjMatrix AprojDat = ProjMatrix(dbin, &mesh);
  ProjMatrix AprojOut = ProjMatrix(dbout, &mesh);

  // Invoke the calculation

  VectorDouble propGlob = dbStatisticsFacies(dbin);
  int ncat = static_cast<int>(propGlob.size());
  OptimCostColored Oc = OptimCostColored(ncat, &Qprop, &AprojDat);

  VectorDouble facies = dbin->getColumnByLocator(ELoc::Z);
  VectorVectorDouble props = Oc.minimize(facies, splits, propGlob, verbose, niter);

  // Loading the resulting results in the output 'dbout'

  int iptr0 = -1;
  VectorDouble propout(dbout->getSampleNumber(true));
  for (int i = 0; i < ncat; i++)
  {
    AprojOut.mesh2point(props[i],propout);
    int iptr = dbout->addColumns(propout,String(),ELoc::UNKNOWN,0,true);
    if (i == 0) iptr0 = iptr;
    namconv.setNamesAndLocators(nullptr, VectorString(), ELoc::UNKNOWN, -1, dbout, iptr,
                                concatenateStrings("-", toString(i + 1)));
  }
  namconv.setLocators(dbout, iptr0, 1, ncat);

  return 0;
}

