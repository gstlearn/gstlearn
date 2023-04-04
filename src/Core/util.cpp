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
#include <Geometry/GeometryHelper.hpp>
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"

#include "LithoRule/Rule.hpp"
#include "Basic/Law.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/PolyLine2D.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceSN.hpp"
#include <string.h>
#include <complex>
#include <cmath>
#include <regex>

/*! \cond */
#define TAB(ix,iy)   (tab[(ix) * ny + (iy)])
#define COORD(i,ip)  (coord[3 * (ip) + (i)])

#define MATTAB(ip,i) (mattab[(ip) * ncolor + (i)])
/*! \endcond */

typedef struct
{
  char keyword[STRING_LENGTH];
  int origin;
  int nrow;
  int ncol;
  double *values;
} Keypair;

typedef struct
{
  int actif;
} Projec_Environ;

typedef struct
{
  int curech;
  int ndim;
  int *nx;
  int *order;
  int *indg;
  double *tab;
} Dim_Loop;

static Projec_Environ PROJEC = { 0 };
static int KEYPAIR_NTAB = 0;
static Keypair *KEYPAIR_TABS = NULL;
static int DISTANCE_NDIM = 0;
static double *DISTANCE_TAB1 = NULL;
static double *DISTANCE_TAB2 = NULL;
static char **LAST_MESSAGE = NULL;
static int NB_LAST_MESSAGE = 0;

/****************************************************************************/
/*!
 **  Toggle the status of the Projection flag
 **
 ** \param[in]  mode Toggle of the projection flag
 ** \li               0   : Switch the flag OFF
 ** \li               1   : Switch the flag ON
 ** \li              -1   : Toggle the flag
 ** \li              else : Do not modify the flag
 **
 *****************************************************************************/
void projec_toggle(int mode)
{
  int projec_actif;

  /* Process the toggling */

  projec_actif = PROJEC.actif;
  if (mode == 1)
    projec_actif = 1;
  else if (mode == 0)
    projec_actif = 0;
  else if (mode == -1) projec_actif = 1 - projec_actif;

  /* Check that the current space is RN */

  bool flag_sphere = (getDefaultSpaceType() == ESpaceType::SN);
  if (flag_sphere && projec_actif)
  {
    messerr("Error when toggling a Projection ON");
    messerr(
          "Definition of a Projection is incompatible with working on a Sphere");
  }
  else
    PROJEC.actif = projec_actif;

  return;
}

/****************************************************************************/
/*!
 **  Returns the projection characteristics
 **
 ** \param[out]  actif activity flag
 **
 *****************************************************************************/
void projec_query(int *actif)

{
  *actif = PROJEC.actif;

  return;
}

/****************************************************************************/
/*!
 **  Print the characteristics of the projection
 **
 *****************************************************************************/
void projec_print(void)

{
  mestitle(1, "Parameters for Projection");
  if (PROJEC.actif)
    message("Projection is switched ON\n");
  else
    message("Projection is switched OFF\n");
  message("Use 'projec.define' to modify previous values\n");
  return;
}

/****************************************************************************/
/*!
 **  Look for an already registered keypair
 **
 ** \return   Rank of the matching item (or -1)
 **
 ** \param[in]  keyword    Keyword
 ** \param[in]  flag_exact 1 if Exact keyword matching is required
 **
 *****************************************************************************/
static int st_match_keypair(const char *keyword, int flag_exact)
{
  Keypair *keypair;
  char keyloc[STRING_LENGTH];

  (void) gslStrcpy(keyloc, keyword);
  string_strip_blanks(keyloc, 0);

  for (int i = 0; i < KEYPAIR_NTAB; i++)
  {
    keypair = &KEYPAIR_TABS[i];
    if (flag_exact)
    {
      if (strcmp(keypair->keyword, keyloc) == 0) return (i);
    }
    else
    {
      if (strstr(keypair->keyword, keyloc) != NULL) return (i);
    }
  }
  return (-1);
}

/****************************************************************************/
/*!
 **  Internal function to find the keypair stack address
 **  or to create a new one if not already existing
 **
 ** \return The address in the stack
 **
 ** \param[in]  keyword Keyword
 **
 ** \remarks If the keypair is new, the arguments 'nrow', 'ncol' and 'origin'
 ** \remarks are set to zero
 ** \remarks Otherwise they are not updated
 **
 *****************************************************************************/
static Keypair* st_get_keypair_address(const char *keyword)

{
  Keypair *keypair;
  char keyloc[STRING_LENGTH];
  int found, flag_new;

  /* Check the length of the keyword */

  if (strlen(keyword) > STRING_LENGTH)
    messageAbort("Keyword %s too long", keyword);

  /* Check if the keyword has already been defined */

  found = st_match_keypair(keyword, 1);
  flag_new = found < 0;

  /* Add a new keypair */

  if (flag_new)
  {
    found = KEYPAIR_NTAB;
    KEYPAIR_NTAB++;
    KEYPAIR_TABS = (Keypair*) realloc((char*) KEYPAIR_TABS,
                                      sizeof(Keypair) * KEYPAIR_NTAB);
  }

  /* Store the attribute (compressing the name and suppressing blanks) */

  keypair = &KEYPAIR_TABS[found];
  (void) gslStrcpy(keyloc, keyword);
  string_strip_blanks(keyloc, 0);
  (void) gslStrcpy(keypair->keyword, keyloc);

  /* Initialize the attributes (for a new keypair) */

  if (flag_new)
  {
    keypair->origin = 0;
    keypair->nrow = 0;
    keypair->ncol = 0;
    keypair->values = NULL;
  }

  return (keypair);
}

/****************************************************************************/
/*!
 **  Internal function to copy or check the attributes (append)
 **
 ** \param[in]  keypair        Keypair structure
 ** \param[in]  mode           0 for creation and 1 for appending
 ** \param[in]  origin         1 from C; 2 from R
 ** \param[in]  ncol           Number of columns
 **
 ** \remarks The arguments 'ncol' and 'origin' are updated.
 ** \remarks Conversely, the argument 'nrow' is not updated here
 **
 *****************************************************************************/
static void st_keypair_attributes(Keypair *keypair,
                                  int mode,
                                  int origin,
                                  int /*nrow*/,
                                  int ncol)
{
  /* Dispatch */

  if (mode == 0)
  {
    // Free the array if attributes are different

    if (keypair->values != NULL)
    {
      if (keypair->ncol != ncol)
      {
        free((char*) keypair->values);
        keypair->values = nullptr;
      }
    }

    // Creation

    keypair->origin = origin;
    keypair->ncol = ncol;
  }
  else
  {

    // Append

    if (keypair->origin == 0 && keypair->ncol == 0)
    {
      keypair->origin = origin;
      keypair->ncol = ncol;
    }
    else
    {
      if (keypair->origin != origin || keypair->ncol != ncol)
        messageAbort(
            "Keypair append cannot change origin or number of columns");
    }
  }
}

/****************************************************************************/
/*!
 **  Internal function to allocate the storage of a keypair
 **
 ** \param[in]  keypair        Keypair structure
 ** \param[in]  nrow           Number of rows
 ** \param[in]  ncol           Number of columns
 **
 *****************************************************************************/
static void st_keypair_allocate(Keypair *keypair, int nrow, int ncol)
{
  int old_size, new_size;

  new_size = nrow * ncol;
  old_size = keypair->nrow * keypair->ncol;

  // If dimensions are unchanged, do nothing 

  if (new_size == old_size && keypair->values != NULL) return;

  // Dimensions are different

  if (old_size == 0)
  {

    // The old dimensions are null, allocate the contents

    keypair->values = (double*) malloc(sizeof(double) * new_size);
  }
  else
  {

    // The old_dimensions are non zero, reallocate the contents

    keypair->values = (double*) realloc((char*) keypair->values,
                                        sizeof(double) * new_size);
  }

  // Ultimate check that allocaiton has been performed correctly

  if (keypair->values == NULL) messageAbort("Keyword allocation failed");

  // Set the number of rows

  keypair->nrow = nrow;
}

/****************************************************************************/
/*!
 **  Internal function to copy the contents of values into he keypair
 **
 ** \param[in]  keypair        Keypair structure
 ** \param[in]  type           1 for integer, 2 for double
 ** \param[in]  start          Staring address within 'values' in keypair
 ** \param[in]  values         Array to be copied
 **
 *****************************************************************************/
static void st_keypair_copy(Keypair *keypair, int type, int start, void *values)
{
  int *icopy, size;
  double *rcopy;

  size = keypair->nrow * keypair->ncol;
  if (type == 1)
  {
    icopy = (int*) values;
    for (int i = 0; i < size; i++)
      keypair->values[i + start] = icopy[i];
  }
  else
  {
    rcopy = (double*) values;
    for (int i = 0; i < size; i++)
      keypair->values[i + start] = rcopy[i];
  }
}

/****************************************************************************/
/*!
 **  Deposit a keypair (double values)
 **
 ** \param[in]  keyword        Keyword
 ** \param[in]  origin         1 from C; 2 from R
 ** \param[in]  nrow           Number of rows
 ** \param[in]  ncol           Number of columns
 ** \param[in]  values         Array of values (Dimension: nrow * ncol)
 **
 ** \remarks All keypair related function use malloc (rather than mem_alloc)
 ** \remarks not to show up in the memory leak calculations
 **
 *****************************************************************************/
void set_keypair(const char *keyword,
                 int origin,
                 int nrow,
                 int ncol,
                 const double *values)
{
  Keypair *keypair;

  /* Get the Keypair address */

  keypair = st_get_keypair_address(keyword);

  /* Store the attributes */

  st_keypair_attributes(keypair, 0, origin, nrow, ncol);

  /* Allocation */

  st_keypair_allocate(keypair, nrow, ncol);

  /* Copy the values */

  st_keypair_copy(keypair, 2, 0, (void*) values);

  return;
}

/****************************************************************************/
/*!
 **  Deposit a keypair (double values) - Append to existing array
 **
 ** \param[in]  keyword        Keyword
 ** \param[in]  origin         1 from C; 2 from R
 ** \param[in]  nrow           Number of rows
 ** \param[in]  ncol           Number of columns
 ** \param[in]  values         Array of values (Dimension: nrow * ncol)
 **
 ** \remarks Appending extends the number of rows but keeps the number of
 ** \remarks columns and the origin unchanged... otherwise fatal error is issued
 ** \remarks All keypair related function use realloc (rather than mem_realloc)
 ** \remarks not to show up in the memory leak calculations
 **
 *****************************************************************************/
void app_keypair(const char *keyword,
                 int origin,
                 int nrow,
                 int ncol,
                 double *values)
{
  Keypair *keypair;
  int start, newrow;

  /* Get the Keypair address */

  keypair = st_get_keypair_address(keyword);

  // If keypair already exists, check that column and origin are unchanged

  if (keypair->ncol > 0)
  {
    if (keypair->ncol != ncol || keypair->origin != origin)
      messageAbort("In 'app_keypair', ncol and origin must be unchaged");
  }

  /* Check the attributes consistency */

  st_keypair_attributes(keypair, 1, origin, nrow, ncol);

  /* Allocation */

  start = ncol * keypair->nrow;
  newrow = nrow + keypair->nrow;
  st_keypair_allocate(keypair, newrow, ncol);

  /* Copy the values */

  st_keypair_copy(keypair, 2, start, (void*) values);
  return;
}

/****************************************************************************/
/*!
 **  Deposit a keypair (integer values)
 **
 ** \param[in]  keyword        Keyword
 ** \param[in]  origin         1 from C; 2 from R
 ** \param[in]  nrow           Number of rows
 ** \param[in]  ncol           Number of columns
 ** \param[in]  values         Array of values (Dimension: nrow * ncol)
 **
 ** \remarks All keypair related function use malloc (rather than mem_alloc)
 ** \remarks not to show up in the memory leak calculations
 **
 *****************************************************************************/
void set_keypair_int(const char *keyword,
                     int origin,
                     int nrow,
                     int ncol,
                     int *values)
{
  Keypair *keypair;

  /* Get the Keypair address */

  keypair = st_get_keypair_address(keyword);

  /* Store the attributes */

  st_keypair_attributes(keypair, 0, origin, nrow, ncol);

  /* Allocation */

  st_keypair_allocate(keypair, nrow, ncol);

  /* Copy the values */

  st_keypair_copy(keypair, 1, 0, (void*) values);
  return;
}

/****************************************************************************/
/*!
 **  Deposit a keypair (doubleinteger values) - Append to existing array
 **
 ** \param[in]  keyword        Keyword
 ** \param[in]  origin         1 from C; 2 from R
 ** \param[in]  nrow           Number of rows
 ** \param[in]  ncol           Number of columns
 ** \param[in]  values         Array of values (Dimension: nrow * ncol)
 **
 ** \remarks Appending extends the number of rows but keeps the number of
 ** \remarks columns and the origin unchanged ... otherwise fatal error is issued
 ** \remarks All keypair related function use realloc (rather than mem_realloc)
 ** \remarks not to show up in the memory leak calculations
 **
 *****************************************************************************/
void app_keypair_int(const char *keyword,
                     int origin,
                     int nrow,
                     int ncol,
                     int *values)
{
  Keypair *keypair;
  int newrow, start;

  /* Get the Keypair address */

  keypair = st_get_keypair_address(keyword);

  // If keypair already exists, check that column and origin are unchanged

  if (keypair->ncol > 0)
  {
    if (keypair->ncol != ncol || keypair->origin != origin)
      messageAbort("In 'app_keypair_int', ncol and origin must be unchaged");
  }

  /* Check the attributes consistency */

  st_keypair_attributes(keypair, 1, origin, nrow, ncol);

  /* Allocation */

  start = ncol * keypair->nrow;
  newrow = nrow + keypair->nrow;
  st_keypair_allocate(keypair, newrow, ncol);

  /* Copy the values */

  st_keypair_copy(keypair, 1, start, (void*) values);
  return;
}

/****************************************************************************/
/*!
 **  Delete a keypair
 **
 ** \param[in]  indice    Index of the Keyword to be deleted
 **
 ** \remarks All keypair related function use malloc (rather than mem_alloc)
 ** \remarks not to show up in the memory leak calculations
 **
 *****************************************************************************/
static void del_keypone(int indice)
{
  Keypair *keypair;

  /* Initializations */

  if (indice < 0 || indice >= KEYPAIR_NTAB) return;

  /* Delete the current keypair */

  keypair = &KEYPAIR_TABS[indice];
  free((char*) keypair->values);
  keypair->values = nullptr;

  /* Shift all subsequent keypairs */

  for (int i = indice + 1; i < KEYPAIR_NTAB; i++)
    KEYPAIR_TABS[i - 1] = KEYPAIR_TABS[i];

  KEYPAIR_NTAB--;
  KEYPAIR_TABS = (Keypair*) realloc((char*) KEYPAIR_TABS,
                                    sizeof(Keypair) * KEYPAIR_NTAB);

  return;
}

/****************************************************************************/
/*!
 **  Delete a keypair
 **
 ** \param[in]  keyword    Keyword to be deleted
 ** \param[in]  flag_exact 1 if Exact keyword matching is required
 **
 ** \remarks All keypair related function use malloc (rather than mem_alloc)
 ** \remarks not to show up in the memory leak calculations
 **
 *****************************************************************************/
void del_keypair(const char *keyword, int flag_exact)
{
  int found;

  if (strlen(keyword) > STRING_LENGTH)
    messageAbort("Keyword %s too long", keyword);

  /* Particular case of the keyword "all" */

  if (!strcmp(keyword, "all"))
  {
    for (int i = KEYPAIR_NTAB - 1; i >= 0; i--)
      del_keypone(i);
  }
  else if (!strcmp(keyword, "allC"))
  {
    for (int i = KEYPAIR_NTAB - 1; i >= 0; i--)
      if (KEYPAIR_TABS[i].origin == 1) del_keypone(i);
  }
  else if (!strcmp(keyword, "allR"))
  {
    for (int i = KEYPAIR_NTAB - 1; i >= 0; i--)
      if (KEYPAIR_TABS[i].origin == 2) del_keypone(i);
  }
  else if (flag_exact)
  {

    /* Delete the keyword with an exact match */

    found = st_match_keypair(keyword, 1);
    if (found < 0) return;

    del_keypone(found);
  }
  else
  {

    /* Delete similar keywords */

    while (1)
    {
      found = st_match_keypair(keyword, 0);
      if (found < 0) return;

      del_keypone(found);
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Inquiry the keypair (for a single value)
 **
 ** \return Returned value
 **
 ** \param[in]  keyword        Keyword
 ** \param[in]  valdef         Factory setting value
 **
 ** \remark  This function will returns systematically the default value
 ** \remark  if the targeted keypair contains more than a single value
 **
 *****************************************************************************/
double get_keypone(const char *keyword, double valdef)
{
  int found;
  double *rtab, retval;
  Keypair *keypair;

  /* Check if the keyword has been defined */

  retval = TEST;
  found = st_match_keypair(keyword, 1);
  if (found >= 0)
  {
    keypair = &KEYPAIR_TABS[found];
    rtab = (double*) keypair->values;
    if (keypair->nrow * keypair->ncol == 1) retval = rtab[0];
  }

  /* Returning argument */

  if (FFFF(retval)) retval = valdef;
  return (retval);
}

/****************************************************************************/
/*!
 **  Inquiry the keypair
 **
 ** \return Error returned code
 **
 ** \param[in]  keyword        Keyword
 **
 ** \param[out] nrow           Number of rows
 ** \param[out] ncol           Number of columns
 ** \param[out] values         Array of values attached to the keyword
 **
 ** \remark  The returned array must be freed by the calling function
 **
 ** \remarks All keypair related function use malloc (rather than mem_alloc)
 ** \remarks not to show up in the memory leak calculations
 **
 *****************************************************************************/
int get_keypair(const char *keyword, int *nrow, int *ncol, double **values)
{
  int found, size;
  double *valloc;
  Keypair *keypair;

  /* Check if the keyword has been defined */

  found = st_match_keypair(keyword, 1);
  if (found < 0) return (1);

  /* The key has been encountered */

  keypair = &KEYPAIR_TABS[found];
  *nrow = keypair->nrow;
  *ncol = keypair->ncol;
  size = (*nrow) * (*ncol);

  valloc = (double*) malloc(sizeof(double) * size);
  for (int i = 0; i < size; i++)
    valloc[i] = keypair->values[i];
  *values = valloc;

  return (0);
}

/****************************************************************************/
/*!
 **  Inquiry the keypair (integer values)
 **
 ** \return Error returned code
 **
 ** \param[in]  keyword        Keyword
 **
 ** \param[out] nrow           Number of rows
 ** \param[out] ncol           Number of columns
 ** \param[out] values         Array of values attached to the keyword
 **
 ** \remark  The returned array must be freed by the calling function
 **
 ** \remarks All keypair related function use malloc (rather than mem_alloc)
 ** \remarks not to show up in the memory leak calculations
 **
 *****************************************************************************/
int get_keypair_int(const char *keyword, int *nrow, int *ncol, int **values)
{
  int *valloc, found, size;
  Keypair *keypair;

  /* Check if the keyword has been defined */

  found = st_match_keypair(keyword, 1);
  if (found < 0) return (1);

  /* The key has been encountered */

  keypair = &KEYPAIR_TABS[found];
  *nrow = keypair->nrow;
  *ncol = keypair->ncol;
  size = (*nrow) * (*ncol);

  valloc = (int*) malloc(sizeof(int) * size);
  for (int i = 0; i < size; i++)
    valloc[i] = (int) keypair->values[i];
  *values = valloc;

  return (0);
}

/****************************************************************************/
/*!
 **  Print the list of keypairs
 **
 ** \param[in]  flag_short  1 for a short output
 **
 *****************************************************************************/
void print_keypair(int flag_short)

{
  int i;
  Keypair *keypair;

  if (KEYPAIR_NTAB <= 0)
    message("No binding keypair is defined\n");
  else
    for (i = 0; i < KEYPAIR_NTAB; i++)
    {
      keypair = &KEYPAIR_TABS[i];
      if (flag_short)
      {
        if (keypair->origin == 1)
          message("C ");
        else
          message("R ");
        message("- %s (%d x %d)\n", keypair->keyword, keypair->nrow,
                keypair->ncol);
      }
      else
        print_matrix(keypair->keyword, 0, 0, keypair->ncol, keypair->nrow, NULL,
                     keypair->values);
    }
  return;
}


/****************************************************************************/
/*!
 **  Find the roots of a polynomial of order 2: ax^2 + bx + c = 0
 **
 ** \return Number of real solutions
 **
 ** \param[in]  a,b,c     Coefficients of the polynomial
 **
 ** \param[out] x         Array of real solutions (Dimension: 2)
 **
 ** \remarks When the solution is double, the returned number os 1.
 **
 *****************************************************************************/
int solve_P2(double a, double b, double c, double *x)
{
  double delta;

  if (a == 0.)
  {
    if (b == 0.)
      return (0);
    else
    {
      x[0] = -c / b;
      return (1);
    }
  }
  else
  {

    // Calculate the discriminant

    delta = b * b - 4 * a * c;

    if (delta == 0.)
    {
      x[0] = -b / (2. * a);
      return (1);
    }
    else
    {
      x[0] = (-b + sqrt(delta)) / (2. * a);
      x[0] = (-b - sqrt(delta)) / (2. * a);
      return (2);
    }
  }
}

/****************************************************************************/
/*!
 **  Find the roots of a polynomial of order 3: a*x^3 + b*x^2 + c*x + d = 0
 **
 ** \return Number of real solutions
 **
 ** \param[in]  a,b,c,d   Coefficients of the polynomial
 **
 ** \param[out] x         Array of real solutions (Dimension: 3)
 **
 ** \remarks When the solution is double, the returned number os 1.
 **
 *****************************************************************************/
int solve_P3(double a, double b, double c, double d, double *x)
{
  double delta, p, q, ecart, u, v, s1;
  int k;

  if (a == 0.)
    return (solve_P2(b, c, d, x));
  else
  {

    // Transform into equation: x^3 + p*x + q = 0

    ecart = -b / (3. * a);
    p = -b * b / (3. * a * a) + c / a;
    q = b / (27. * a) * (2. * b * b / (a * a) - 9. * c / a) + d / a;

    // Cardan formula

    delta = -(4. * p * p * p + 27. * q * q);
    if (delta < 0)
    {
      s1 = sqrt(-delta / 27.);
      u = (-q + s1) / 2.;
      u = (u > 0.) ? pow(u, 1. / 3.) : -pow(-u, 1. / 3.);
      v = (-q - s1) / 2.;
      v = (v > 0.) ? pow(v, 1. / 3.) : -pow(-v, 1. / 3.);
      x[0] = ecart + u + v;
      return (1);
    }
    else if (delta == 0.)
    {
      x[0] = ecart + 3. * q / p;
      x[1] = ecart - 3. * q / (2. * p);
      return (2);
    }
    else
    {
      s1 = -(q / 2.) * sqrt(27. / -(p * p * p));
      for (k = 0; k < 3; k++)
        x[k] = ecart + 2. * sqrt(-p / 3.) * cos((acos(s1) + 2. * k * GV_PI) / 3.);
      return (3);
    }
  }
}

/****************************************************************************/
/*!
 **  Calculate the distance between two endpoints
 **
 ** \return Distance value (or TEST if a coordinate is not defined)
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  tab1   Array corresponding to the first endpoint
 ** \param[in]  tab2   Array corresponding to the second endpoint
 **
 *****************************************************************************/
double ut_distance(int ndim, const double *tab1, const double *tab2)
{
  double distance, v1, v2, delta;

  distance = 0.;
  bool flag_sphere = isDefaultSpaceSphere();

  if (flag_sphere)
  {
    /* Case of the spherical coordinates */
    /* Longitude = 1st coord; Latitude = 2nd coord (in degrees) */

    const ASpace* space = getDefaultSpace();
    const SpaceSN* spaceSn = dynamic_cast<const SpaceSN*>(space);
    if (space == nullptr) return TEST;
    double R = spaceSn->getRadius();
    distance = GH::geodeticAngularDistance(tab1[0], tab1[1], tab2[0], tab2[1], R);
  }
  else
  {
    /* Case of the euclidean coordinates */

    for (int idim = 0; idim < ndim; idim++)
    {
      v1 = tab1[idim];
      if (FFFF(v1)) return (TEST);
      v2 = tab2[idim];
      if (FFFF(v2)) return (TEST);
      delta = v1 - v2;
      distance += delta * delta;
    }
    distance = sqrt(distance);
  }
  return (distance);
}

/*****************************************************************************/
/*!
 **  Allocate the necessary arrays for calculating distances
 **  using already allocated arrays
 **
 ** \param[in]  ndim   Space dimension
 **
 ** \param[out] tab1   Array for coordinates of first sample
 ** \param[out] tab2   Array for coordinates of second sample
 **
 ** \remarks This function uses realloc (rather than mem_realloc)
 ** \remarks not to show up in the memory leak calculations
 **
 *****************************************************************************/
void ut_distance_allocated(int ndim, double **tab1, double **tab2)
{
  if (DISTANCE_NDIM < ndim)
  {
    DISTANCE_TAB1 = (double*) realloc((char*) DISTANCE_TAB1,
                                      sizeof(double) * ndim);
    DISTANCE_TAB2 = (double*) realloc((char*) DISTANCE_TAB2,
                                      sizeof(double) * ndim);
    DISTANCE_NDIM = ndim;
  }
  *tab1 = DISTANCE_TAB1;
  *tab2 = DISTANCE_TAB2;
  return;
}


/****************************************************************************/
/*!
 ** Create a VectorDouble for storing an array of double
 **
 ** \return The VectorDouble
 **
 ** \param[in]  ntab      Number of samples
 ** \param[in]  rtab      Array of double values to be loaded
 **
 *****************************************************************************/
VectorDouble util_set_array_double(int ntab, const double *rtab)
{
  if (OptDbg::query(EDbg::INTERFACE)) message("util_set_array_double\n");
  if (ntab <= 0 || rtab == nullptr) return VectorDouble();
  VectorDouble rettab(ntab);
  if (rettab.empty()) return rettab;

  for (int i = 0; i < ntab; i++)
    rettab[i] = rtab[i];

  return rettab;
}

/****************************************************************************/
/*!
 ** Create a VectorInt for storing an array of integer
 **
 ** \return  The VectorInt
 **
 ** \param[in]  ntab      Number of samples
 ** \param[in]  itab      Array of integer values to be loaded
 **
 *****************************************************************************/
VectorInt util_set_array_integer(int ntab, const int *itab)
{
  if (OptDbg::query(EDbg::INTERFACE)) message("util_set_array_integer\n");
  VectorInt rettab(ntab);
  if (ntab <= 0 || itab == nullptr) return rettab;
  for (int i = 0; i < ntab; i++)
    rettab[i] = itab[i];
  return rettab;
}

/****************************************************************************/
/*!
 ** Create a VectorString for storing an array of chars
 **
 ** \return  The VectorString
 **
 ** \param[in]  ntab      Number of samples
 ** \param[in]  names     Array of character values to be loaded
 **
 *****************************************************************************/
VectorString util_set_array_char(int ntab, char **names)
{
  if (OptDbg::query(EDbg::INTERFACE)) message("util_set_array_char\n");
  VectorString rettab(ntab);
  if (names == nullptr) return rettab;
  for (int i = 0; i < ntab; i++)
    rettab[i] = names[i];
  return rettab;
}

/****************************************************************************/
/*!
 **  Deposit a last message
 **
 ** \param[in]  mode           Type of operation
 **                            0 to empty the array of messages
 **                            1 to add the string to the array of messages
 **                           -1 to concatenate the string to the last message
 ** \param[in]  string         Current string
 **
 ** \remarks All keypair related function use malloc (rather than mem_alloc)
 ** \remarks not to show up in the memory leak calculations
 **
 *****************************************************************************/
void set_last_message(int mode, const char *string)
{
  char *address;
  int size, sizaux;

  /* Dispatch */

  switch (mode)
  {
    case 0:
      if (NB_LAST_MESSAGE <= 0) return;
      for (int i = 0; i < NB_LAST_MESSAGE; i++)
      {
        free((char*) LAST_MESSAGE[i]);
        LAST_MESSAGE[i] = nullptr;
      }
      free((char*) LAST_MESSAGE);
      NB_LAST_MESSAGE = 0;
      break;

    case 1:                       // Add string to array of messages
      size = static_cast<int>(strlen(string));
      if (size <= 0) return;

      if (NB_LAST_MESSAGE <= 0)
        LAST_MESSAGE = (char**) malloc(sizeof(char*) * 1);
      else
        LAST_MESSAGE = (char**) realloc((char*) LAST_MESSAGE,
                                        sizeof(char*) * (NB_LAST_MESSAGE + 1));
      LAST_MESSAGE[NB_LAST_MESSAGE] = address = (char*) malloc(size + 1);
      (void) gslStrcpy(address, string);
      address[size] = '\0';
      NB_LAST_MESSAGE++;
      break;

    case -1:                    // Concatenate
      size = static_cast<int>(strlen(string));
      if (size <= 0) return;

      if (NB_LAST_MESSAGE <= 0)
      {
        set_last_message(1, string);
        return;
      }

      sizaux = static_cast<int>(strlen(LAST_MESSAGE[NB_LAST_MESSAGE - 1]));
      LAST_MESSAGE[NB_LAST_MESSAGE - 1] = address = (char*) realloc(
          (char*) LAST_MESSAGE[NB_LAST_MESSAGE - 1], size + sizaux + 2);
      address[sizaux] = ' ';
      (void) gslStrcpy(&address[sizaux + 1], string);
      address[size + sizaux + 1] = '\0';
      break;
  }
}

/****************************************************************************/
/*!
 **  Print the array of last messages
 **
 *****************************************************************************/
void print_last_message(void)
{
  if (NB_LAST_MESSAGE <= 0) return;

  mestitle(0, "Last Message");
  for (int i = 0; i < NB_LAST_MESSAGE; i++)
  {
    message(">>> %s\n", LAST_MESSAGE[i]);
  }
  message("\n");
}

/****************************************************************************/
/*!
 **  Return all the ways to split ncolor into two non-empty subsets
 **
 ** \return Return an array of possibilities
 **
 ** \param[in]  ncolor    Number of colors
 ** \param[in]  flag_half 1 if only half of possibilities must be envisaged
 ** \param[in]  verbose   1 for a verbose option
 **
 ** \param[out] nposs  Number of possibilities
 **
 ** \remarks The calling function must free the returned array.
 ** \remarks The array has 'ncolor' columns and 'ncomb' subsets
 ** \remarks The elements of each row are set to 0 or 1 (subset rank)
 **
 *****************************************************************************/
int* ut_split_into_two(int ncolor, int flag_half, int verbose, int *nposs)
{
  int p, nmax, ncomb, np, lec;
  int *mattab, *comb;

  /* Initializations */

  p = (flag_half) ? static_cast<int>(floor((double) ncolor / 2.)) :
                    ncolor - 1;
  nmax = static_cast<int>(pow(2, ncolor));
  mattab = comb = nullptr;
  np = 0;

  /* Core allocation */

  mattab = (int*) mem_alloc(sizeof(int) * ncolor * nmax, 1);
  for (int i = 0; i < ncolor * nmax; i++)
    mattab[i] = 0;

  for (int nsub = 1; nsub <= p; nsub++)
  {
    comb = ut_combinations(ncolor, nsub, &ncomb);
    lec = 0;
    for (int i = 0; i < ncomb; i++)
    {
      for (int j = 0; j < nsub; j++, lec++)
        MATTAB(np,comb[lec]-1) = 1;
      np++;
    }
  }
  comb = (int*) mem_free((char* ) comb);

  /* Resize */

  mattab = (int*) mem_realloc((char* ) mattab, sizeof(int) * ncolor * np, 1);
  *nposs = np;

  /* Verbose option */

  if (verbose)
  {
    message("Initial number of values = %d (Half=%d)\n", ncolor, flag_half);
    lec = 0;
    for (int i = 0; i < np; i++)
    {
      for (int j = 0; j < ncolor; j++, lec++)
        message(" %d", mattab[lec]);
      message("\n");
    }
  }
  return (mattab);
}

/****************************************************************************/
/*!
 **   Convert std::string into a char *
 **
 ** \return  Pointer to the returned array of characters
 **
 ** \param[in]  s        Input VectorString
 **
 *****************************************************************************/
char* convert(const std::string &s)
{
  char *pc = new char[s.size() + 1];
  (void) gslStrcpy(pc, s.c_str());
  return pc;
}

/****************************************************************************/
/*!
 **   Convert VectorString into a std::vector<char *> structure
 **
 ** \return  Pointer to the returned array of characters
 **
 ** \param[in]  vs        Input VectorString
 **
 *****************************************************************************/
std::vector<char*> util_vs_to_vs(VectorString vs)
{
  std::vector<char*> vc;
  std::transform(vs.begin(), vs.end(), std::back_inserter(vc), convert);
  return vc;
}
