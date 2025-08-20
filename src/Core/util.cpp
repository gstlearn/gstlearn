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
#include "Geometry/GeometryHelper.hpp"
#include "Space/ASpace.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpaceSN.hpp"
#include "geoslib_old_f.h"

#include <cmath>
#include <cstring>

/*! \cond */
#define TAB(ix, iy)  (tab[(ix) * ny + (iy)])
#define COORD(i, ip) (coord[3 * (ip) + (i)])

#define MATTAB(ip, i) (mattab[(ip) * ncolor + (i)])
/*! \endcond */

namespace gstlrn
{
typedef struct
{
  char keyword[STRING_LENGTH];
  Id origin;
  Id nrow;
  Id ncol;
  VectorDouble values;
} Keypair;

static Id KEYPAIR_NTAB = 0;
static std::vector<Keypair> KEYPAIR_TABS;
static Id DISTANCE_NDIM = 0;
static VectorDouble DISTANCE_TAB1;
static VectorDouble DISTANCE_TAB2;

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
static Id st_match_keypair(const char* keyword, Id flag_exact)
{
  Keypair* keypair;
  char keyloc[STRING_LENGTH];

  (void)gslStrcpy(keyloc, STRING_LENGTH, keyword);
  string_strip_blanks(keyloc, 0);

  for (Id i = 0; i < KEYPAIR_NTAB; i++)
  {
    keypair = &KEYPAIR_TABS[i];
    if (flag_exact)
    {
      if (strcmp(keypair->keyword, keyloc) == 0) return (i);
    }
    else
    {
      if (strstr(keypair->keyword, keyloc) != nullptr) return (i);
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
static Keypair* st_get_keypair_address(const char* keyword)

{
  Keypair* keypair;
  char keyloc[STRING_LENGTH];
  Id found, flag_new;

  /* Check the length of the keyword */

  if (strlen(keyword) > STRING_LENGTH)
    messageAbort("Keyword %s too long", keyword);

  /* Check if the keyword has already been defined */

  found    = st_match_keypair(keyword, 1);
  flag_new = found < 0;

  /* Add a new keypair */

  if (flag_new)
  {
    found = KEYPAIR_NTAB;
    KEYPAIR_NTAB++;
    KEYPAIR_TABS.resize(KEYPAIR_NTAB);
  }

  /* Store the attribute (compressing the name and suppressing blanks) */

  keypair = &KEYPAIR_TABS[found];
  (void)gslStrcpy(keyloc, STRING_LENGTH, keyword);
  string_strip_blanks(keyloc, 0);
  (void)gslStrcpy(keypair->keyword, STRING_LENGTH, keyloc);

  /* Initialize the attributes (for a new keypair) */

  if (flag_new)
  {
    keypair->origin = 0;
    keypair->nrow   = 0;
    keypair->ncol   = 0;
    keypair->values.clear();
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
static void st_keypair_attributes(Keypair* keypair,
                                  Id mode,
                                  Id origin,
                                  Id /*nrow*/,
                                  Id ncol)
{
  /* Dispatch */

  if (mode == 0)
  {
    // Free the array if attributes are different

    if (!keypair->values.empty())
    {
      if (keypair->ncol != ncol)
      {
        keypair->values.clear();
      }
    }

    // Creation

    keypair->origin = origin;
    keypair->ncol   = ncol;
  }
  else
  {

    // Append

    if (keypair->origin == 0 && keypair->ncol == 0)
    {
      keypair->origin = origin;
      keypair->ncol   = ncol;
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
static void st_keypair_allocate(Keypair* keypair, Id nrow, Id ncol)
{

  const auto new_size = nrow * ncol;
  keypair->values.resize(new_size);

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
static void st_keypair_copy(Keypair* keypair, Id type, Id start, void* values)
{
  Id *icopy, size;
  double* rcopy;

  size = keypair->nrow * keypair->ncol;
  if (type == 1)
  {
    icopy = static_cast<Id*>(values);
    for (Id i = 0; i < size; i++)
      keypair->values[i + start] = icopy[i];
  }
  else
  {
    rcopy = static_cast<double*>(values);
    for (Id i = 0; i < size; i++)
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
 *****************************************************************************/
void set_keypair(const char* keyword,
                 Id origin,
                 Id nrow,
                 Id ncol,
                 const double* values)
{
  Keypair* keypair;

  /* Get the Keypair address */

  keypair = st_get_keypair_address(keyword);

  /* Store the attributes */

  st_keypair_attributes(keypair, 0, origin, nrow, ncol);

  /* Allocation */

  st_keypair_allocate(keypair, nrow, ncol);

  /* Copy the values */

  st_keypair_copy(keypair, 2, 0, static_cast<void*>(const_cast<double*>(values)));
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
 **
 *****************************************************************************/
void app_keypair(const char* keyword,
                 Id origin,
                 Id nrow,
                 Id ncol,
                 double* values)
{
  Keypair* keypair;
  Id start, newrow;

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

  start  = ncol * keypair->nrow;
  newrow = nrow + keypair->nrow;
  st_keypair_allocate(keypair, newrow, ncol);

  /* Copy the values */

  st_keypair_copy(keypair, 2, start, static_cast<void*>(values));
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
 *****************************************************************************/
void set_keypair_int(const char* keyword,
                     Id origin,
                     Id nrow,
                     Id ncol,
                     Id* values)
{
  Keypair* keypair;

  /* Get the Keypair address */

  keypair = st_get_keypair_address(keyword);

  /* Store the attributes */

  st_keypair_attributes(keypair, 0, origin, nrow, ncol);

  /* Allocation */

  st_keypair_allocate(keypair, nrow, ncol);

  /* Copy the values */

  st_keypair_copy(keypair, 1, 0, static_cast<void*>(values));
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
 **
 *****************************************************************************/
void app_keypair_int(const char* keyword,
                     Id origin,
                     Id nrow,
                     Id ncol,
                     Id* values)
{
  Keypair* keypair;
  Id newrow, start;

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

  start  = ncol * keypair->nrow;
  newrow = nrow + keypair->nrow;
  st_keypair_allocate(keypair, newrow, ncol);

  /* Copy the values */

  st_keypair_copy(keypair, 1, start, static_cast<void*>(values));
}

/****************************************************************************/
/*!
 **  Delete a keypair
 **
 ** \param[in]  indice    Index of the Keyword to be deleted
 **
 *****************************************************************************/
static void del_keypone(Id indice)
{
  Keypair* keypair;

  /* Initializations */

  if (indice < 0 || indice >= KEYPAIR_NTAB) return;

  /* Delete the current keypair */

  keypair = &KEYPAIR_TABS[indice];
  keypair->values.clear();

  /* Shift all subsequent keypairs */

  for (Id i = indice + 1; i < KEYPAIR_NTAB; i++)
    KEYPAIR_TABS[i - 1] = KEYPAIR_TABS[i];

  KEYPAIR_NTAB--;
  KEYPAIR_TABS.resize(KEYPAIR_NTAB);
}

/****************************************************************************/
/*!
 **  Delete a keypair
 **
 ** \param[in]  keyword    Keyword to be deleted
 ** \param[in]  flag_exact 1 if Exact keyword matching is required
 **
 *****************************************************************************/
void del_keypair(const char* keyword, Id flag_exact)
{
  Id found;

  if (strlen(keyword) > STRING_LENGTH)
    messageAbort("Keyword %s too long", keyword);

  /* Particular case of the keyword "all" */

  if (!strcmp(keyword, "all"))
  {
    for (Id i = KEYPAIR_NTAB - 1; i >= 0; i--)
      del_keypone(i);
  }
  else if (!strcmp(keyword, "allC"))
  {
    for (Id i = KEYPAIR_NTAB - 1; i >= 0; i--)
      if (KEYPAIR_TABS[i].origin == 1) del_keypone(i);
  }
  else if (!strcmp(keyword, "allR"))
  {
    for (Id i = KEYPAIR_NTAB - 1; i >= 0; i--)
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
double get_keypone(const char* keyword, double valdef)
{
  Id found;
  double retval;
  Keypair* keypair;

  /* Check if the keyword has been defined */

  retval = TEST;
  found  = st_match_keypair(keyword, 1);
  if (found >= 0)
  {
    keypair          = &KEYPAIR_TABS[found];
    const auto& rtab = keypair->values;
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
 *****************************************************************************/
Id get_keypair(const char* keyword, Id* nrow, Id* ncol, VectorDouble& values)
{
  Id found, size;
  Keypair* keypair;

  /* Check if the keyword has been defined */

  found = st_match_keypair(keyword, 1);
  if (found < 0) return (1);

  /* The key has been encountered */

  keypair = &KEYPAIR_TABS[found];
  *nrow   = keypair->nrow;
  *ncol   = keypair->ncol;
  size    = (*nrow) * (*ncol);

  values.resize(size);
  for (Id i = 0; i < size; i++)
    values[i] = keypair->values[i];

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
 *****************************************************************************/
Id get_keypair_int(const char* keyword, Id* nrow, Id* ncol, VectorInt& values)
{
  Id found, size;
  Keypair* keypair;

  /* Check if the keyword has been defined */

  found = st_match_keypair(keyword, 1);
  if (found < 0) return (1);

  /* The key has been encountered */

  keypair = &KEYPAIR_TABS[found];
  *nrow   = keypair->nrow;
  *ncol   = keypair->ncol;
  size    = (*nrow) * (*ncol);

  values.resize(size);
  for (Id i = 0; i < size; i++)
    values[i] = static_cast<Id>(keypair->values[i]);

  return (0);
}

/****************************************************************************/
/*!
 **  Print the list of keypairs
 **
 ** \param[in]  flag_short  1 for a short output
 **
 *****************************************************************************/
void print_keypair(Id flag_short)

{
  Id i;
  Keypair* keypair;

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
                     keypair->values.data());
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
double ut_distance(Id ndim, const double* tab1, const double* tab2)
{
  double distance, v1, v2, delta;

  distance         = 0.;
  bool flag_sphere = isDefaultSpaceSphere();

  if (flag_sphere)
  {
    /* Case of the spherical coordinates */
    /* Longitude = 1st coord; Latitude = 2nd coord (in degrees) */

    const ASpace* space = getDefaultSpaceSh().get();
    const auto* spaceSn = dynamic_cast<const SpaceSN*>(space);
    if (space == nullptr) return TEST;
    double R = spaceSn->getRadius();
    distance = GH::geodeticAngularDistance(tab1[0], tab1[1], tab2[0], tab2[1], R);
  }
  else
  {
    /* Case of the euclidean coordinates */

    for (Id idim = 0; idim < ndim; idim++)
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
 *****************************************************************************/
void ut_distance_allocated(Id ndim, double** tab1, double** tab2)
{
  if (DISTANCE_NDIM < ndim)
  {
    DISTANCE_TAB1.resize(ndim);
    DISTANCE_TAB2.resize(ndim);
    DISTANCE_NDIM = ndim;
  }
  *tab1 = DISTANCE_TAB1.data();
  *tab2 = DISTANCE_TAB2.data();
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
VectorInt ut_split_into_two(Id ncolor, Id flag_half, Id verbose, Id* nposs)
{
  Id p, nmax, ncomb, np, lec;
  VectorInt mattab;

  /* Initializations */

  p    = (flag_half) ? static_cast<Id>(floor(static_cast<double>(ncolor) / 2.)) : ncolor - 1;
  nmax = static_cast<Id>(pow(2, ncolor));
  np   = 0;

  /* Core allocation */

  mattab.resize(ncolor * nmax);

  for (Id nsub = 1; nsub <= p; nsub++)
  {
    VectorInt comb = ut_combinations(ncolor, nsub, &ncomb);
    lec            = 0;
    for (Id i = 0; i < ncomb; i++)
    {
      for (Id j = 0; j < nsub; j++, lec++)
        MATTAB(np, comb[lec] - 1) = 1;
      np++;
    }
  }

  /* Resize */

  mattab.resize(ncolor * np);
  *nposs = np;

  /* Verbose option */

  if (verbose)
  {
    message("Initial number of values = %d (Half=%d)\n", ncolor, flag_half);
    lec = 0;
    for (Id i = 0; i < np; i++)
    {
      for (Id j = 0; j < ncolor; j++, lec++)
        message(" %d", mattab[lec]);
      message("\n");
    }
  }
  return (mattab);
}

} // namespace gstlrn
