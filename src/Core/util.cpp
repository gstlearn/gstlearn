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
#include "Basic/Law.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "geoslib_e.h"
#include "geoslib_old_f.h"
//#include <tr1/cmath>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex>
#include <cmath>
#include <regex>

/*! \cond */
#define ROT(i,j)     (rot[(i) * ndim + (j)])
#define TAB(ix,iy)   (tab[(ix) * ny + (iy)])
#define COORD(i,ip)  (coord[3 * (ip) + (i)])
#define RCOORD(i,ip) (R_coor->coor[3 * (ip) + (i)])
#define MATTAB(ip,i) (mattab[(ip) * ncolor + (i)])
#define DBG_NUMBER  15
#define MINI        10
#define LSTACK    1000
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
  int status;
  char keyword[20];
  char comment[STRING_LENGTH];
} Debug;

static Debug DBG[DBG_NUMBER] = { { 0, "interface", "Communication with interface" },
                                 { 0, "db",        "Data Base Management" },
                                 { 0, "nbgh",      "Neighborhood Management" },
                                 { 0, "model",     "Model Management" },
                                 { 0, "kriging",   "Kriging Operations" },
                                 { 0, "simulate",  "Simulations" },
                                 { 0, "results",   "Kriging Results" },
                                 { 0, "variogram", "Variogram calculations" },
                                 { 0, "converge",  "Convergence test" },
                                 { 0, "condexp",   "Conditional Expectation" },
                                 { 0, "bayes",     "Bayesian Estimation" },
                                 { 0, "morpho",    "Morphological Operations" },
                                 { 0, "props",     "Proportions or Intensities" },
                                 { 0, "upscale",   "Upscaling" },
                                 { 0, "spde",      "S.P.D.E" } };

typedef struct
{
  int actif;
} Projec_Environ;

typedef struct
{
  int index;
  int reference;
} Debug_Environ;

typedef struct
{
  int flag_sphere;
  double radius;
} Variety_Environ;

typedef struct
{
  int ntri;
  double *coor;
} Reg_Coor;

typedef struct
{
  int curech;
  int ndim;
  int *nx;
  int *order;
  int *indg;
  double *tab;
} Dim_Loop;

static double (*LEGENDRE_SPHPLM)(int, int, double) = NULL;
static double (*LEGENDRE_PL)(int, double) = NULL;

static Projec_Environ PROJEC = { 0 };
static Debug_Environ DBGENV = { 0, 0 };
static Variety_Environ VARIETY = { 0, 0. };
static int KEYPAIR_NTAB = 0;
static Keypair *KEYPAIR_TABS = NULL;
static char QUESTION[STRING_LENGTH];
static int DISTANCE_NDIM = 0;
static double *DISTANCE_TAB1 = NULL;
static double *DISTANCE_TAB2 = NULL;
static char INSTR1[STRING_LENGTH];
static char INSTR2[STRING_LENGTH];
static char **LAST_MESSAGE = NULL;
static int NB_LAST_MESSAGE = 0;


/*****************************************************************************/
/*!
 **  Returns the unique occurrence of values in a vector of values
 **
 ** \param[in]  ntab   Number of values
 ** \param[in]  tab    Array of values
 **
 ** \param[out] neff   Actual number of different values
 **
 ** \remark  The 'neff' values are placed at the beginning of 'tab' in output
 **
 *****************************************************************************/
GEOSLIB_API void ut_tab_unique(int ntab, double *tab, int *neff)
{
  int ecr;
  double value;

  /* Sorting the values */

  ut_sort_double(0, ntab, NULL, tab);

  /* Search for unique occurences */

  ecr = 1;
  value = tab[0];
  for (int lec = 1; lec < ntab; lec++)
  {
    if (tab[lec] == value) continue;
    value = tab[lec];
    tab[ecr++] = value;
  }
  *neff = ecr;
}

/*****************************************************************************/
/*!
 **  Sorts the (double) array value() and the array ind()
 **  in the ascending order of value
 **
 ** \param[in]  safe   1 if the value array if preserved
 **                    0 if the value array is also sorted
 ** \param[in]  nech   number of samples
 **
 ** \param[out] ind    output int array
 ** \param[out] value  input and output array
 **
 ** \remark  If ind = NULL, ind is ignored
 **
 *****************************************************************************/
GEOSLIB_API void ut_sort_double(int safe, int nech, int *ind, double *value)
{
  static int LISTE_L[LSTACK];
  static int LISTE_R[LSTACK];
  int i, j, p, l, r, pstack, inddev, inddeu;
  double *tab, tablev, tableu;

  /* Initialization */

  inddev = inddeu = 0;
  if (safe)
  {
    tab = (double *) mem_alloc(sizeof(double) * nech, 1);
    for (i = 0; i < nech; i++)
      tab[i] = value[i];
  }
  else
    tab = value;

  /* Segmentation */

  if (nech > MINI)
  {
    if (ind)
    {
      pstack = 0;
      LISTE_L[pstack] = 0;
      LISTE_R[pstack] = nech - 1;

      while (pstack >= 0)
      {
        l = LISTE_L[pstack];
        r = LISTE_R[pstack];
        while (r - l + 1 > MINI)
        {
          i = l;
          j = r + 1;
          p = (int) ((double) (l + r) / 2.);

          if (tab[p] < tab[l])
          {
            tablev = tab[p];
            tab[p] = tab[l];
            tab[l] = tablev;
            inddev = ind[p];
            ind[p] = ind[l];
            ind[l] = inddev;
          }
          if (tab[l] < tab[r])
          {
            if (tab[p] > tab[r])
            {
              tablev = tab[r];
              tab[r] = tab[l];
              tab[l] = tablev;
              inddev = ind[r];
              ind[r] = ind[l];
              ind[l] = inddev;
            }
            else
            {
              tablev = tab[p];
              tab[p] = tab[l];
              tab[l] = tablev;
              inddev = ind[p];
              ind[p] = ind[l];
              ind[l] = inddev;
            }
          }
          else
          {
            tablev = tab[l];
            inddev = ind[l];
          }

          while (1)
          {
            i++;
            while (tab[i] < tablev)
              i++;
            j--;
            while (tab[j] > tablev)
              j--;
            if (j <= i) break;
            tableu = tab[i];
            tab[i] = tab[j];
            tab[j] = tableu;
            inddeu = ind[i];
            ind[i] = ind[j];
            ind[j] = inddeu;
          }
          tab[l] = tab[j];
          tab[j] = tablev;
          ind[l] = ind[j];
          ind[j] = inddev;

          if (r - j > j - l)
          {
            if (j - l > MINI)
            {
              LISTE_L[pstack] = l;
              LISTE_R[pstack] = j - 1;
              pstack++;
              if (pstack >= LSTACK) messageAbort("Stack Overflow");
            }
            l = j + 1;
          }
          else
          {
            if (r - j > MINI)
            {
              LISTE_L[pstack] = j + 1;
              LISTE_R[pstack] = r;
              pstack++;
              if (pstack >= LSTACK) messageAbort("Stack Overflow");
            }
            r = j - 1;
          }
        }
        pstack--;
      }
    }
    else
    {
      pstack = 0;
      LISTE_L[pstack] = 0;
      LISTE_R[pstack] = nech - 1;

      while (pstack >= 0)
      {
        l = LISTE_L[pstack];
        r = LISTE_R[pstack];
        while (r - l + 1 > MINI)
        {
          i = l;
          j = r + 1;
          p = (int) ((double) (l + r) / 2.);

          if (tab[p] < tab[l])
          {
            tablev = tab[p];
            tab[p] = tab[l];
            tab[l] = tablev;
          }
          if (tab[l] < tab[r])
          {
            if (tab[p] > tab[r])
            {
              tablev = tab[r];
              tab[r] = tab[l];
              tab[l] = tablev;
            }
            else
            {
              tablev = tab[p];
              tab[p] = tab[l];
              tab[l] = tablev;
            }
          }
          else
          {
            tablev = tab[l];
          }

          while (1)
          {
            i++;
            while (tab[i] < tablev)
              i++;
            j--;
            while (tab[j] > tablev)
              j--;
            if (j <= i) break;
            tableu = tab[i];
            tab[i] = tab[j];
            tab[j] = tableu;
          }
          tab[l] = tab[j];
          tab[j] = tablev;

          if (r - j > j - l)
          {
            if (j - l > MINI)
            {
              LISTE_L[pstack] = l;
              LISTE_R[pstack] = j - 1;
              pstack++;
              if (pstack >= LSTACK) messageAbort("Stack Overflow");
            }
            l = j + 1;
          }
          else
          {
            if (r - j > MINI)
            {
              LISTE_L[pstack] = j + 1;
              LISTE_R[pstack] = r;
              pstack++;
              if (pstack >= LSTACK) messageAbort("Stack Overflow");
            }
            r = j - 1;
          }
        }
        pstack--;
      }
    }
  }

  /* Final sorting */

  for (j = 1; j < nech; j++)
  {
    if (tab[j - 1] > tab[j])
    {
      tablev = tab[j];
      if (ind) inddev = ind[j];
      i = j;
      while (i > 0)
      {
        if (tab[i - 1] <= tablev) break;
        tab[i] = tab[i - 1];
        if (ind) ind[i] = ind[i - 1];
        i--;
      }
      tab[i] = tablev;
      if (ind) ind[i] = inddev;
    }
  }

  if (safe) tab = (double *) mem_free((char * ) tab);
  return;
}

/*****************************************************************************/
/*!
 **  Sorts the (integer) array value() and the array ind()
 **  in the ascending order of value
 **
 ** \param[in]  safe  1 if the value array if preserved
 **                   0 if the value array is also sorted
 ** \param[in]  nech  number of samples
 **
 ** \param[out] ind    output int array
 ** \param[out] value  input and output array
 **
 ** \remark  If ind = NULL, ind is ignored
 **
 *****************************************************************************/
GEOSLIB_API void ut_sort_int(int safe, int nech, int *ind, int *value)
{
  static int LISTE_L[LSTACK];
  static int LISTE_R[LSTACK];
  int i, j, p, l, r, pstack, inddev, inddeu;
  int *tab, tablev, tableu;

  /* Initialization */

  inddev = inddeu = 0;
  if (safe)
  {
    tab = (int *) mem_alloc(sizeof(int) * nech, 1);
    for (i = 0; i < nech; i++)
      tab[i] = value[i];
  }
  else
    tab = value;

  /* Segmentation */

  if (nech > MINI)
  {
    if (ind)
    {
      pstack = 0;
      LISTE_L[pstack] = 0;
      LISTE_R[pstack] = nech - 1;

      while (pstack >= 0)
      {
        l = LISTE_L[pstack];
        r = LISTE_R[pstack];
        while (r - l + 1 > MINI)
        {
          i = l;
          j = r + 1;
          p = (int) ((double) (l + r) / 2.);

          if (tab[p] < tab[l])
          {
            tablev = tab[p];
            tab[p] = tab[l];
            tab[l] = tablev;
            inddev = ind[p];
            ind[p] = ind[l];
            ind[l] = inddev;
          }
          if (tab[l] < tab[r])
          {
            if (tab[p] > tab[r])
            {
              tablev = tab[r];
              tab[r] = tab[l];
              tab[l] = tablev;
              inddev = ind[r];
              ind[r] = ind[l];
              ind[l] = inddev;
            }
            else
            {
              tablev = tab[p];
              tab[p] = tab[l];
              tab[l] = tablev;
              inddev = ind[p];
              ind[p] = ind[l];
              ind[l] = inddev;
            }
          }
          else
          {
            tablev = tab[l];
            inddev = ind[l];
          }

          while (1)
          {
            i++;
            while (tab[i] < tablev)
              i++;
            j--;
            while (tab[j] > tablev)
              j--;
            if (j <= i) break;
            tableu = tab[i];
            tab[i] = tab[j];
            tab[j] = tableu;
            inddeu = ind[i];
            ind[i] = ind[j];
            ind[j] = inddeu;
          }
          tab[l] = tab[j];
          tab[j] = tablev;
          ind[l] = ind[j];
          ind[j] = inddev;

          if (r - j > j - l)
          {
            if (j - l > MINI)
            {
              LISTE_L[pstack] = l;
              LISTE_R[pstack] = j - 1;
              pstack++;
              if (pstack >= LSTACK) messageAbort("Stack Overflow");
            }
            l = j + 1;
          }
          else
          {
            if (r - j > MINI)
            {
              LISTE_L[pstack] = j + 1;
              LISTE_R[pstack] = r;
              pstack++;
              if (pstack >= LSTACK) messageAbort("Stack Overflow");
            }
            r = j - 1;
          }
        }
        pstack--;
      }
    }
    else
    {
      pstack = 0;
      LISTE_L[pstack] = 0;
      LISTE_R[pstack] = nech - 1;

      while (pstack >= 0)
      {
        l = LISTE_L[pstack];
        r = LISTE_R[pstack];
        while (r - l + 1 > MINI)
        {
          i = l;
          j = r + 1;
          p = (int) ((double) (l + r) / 2.);

          if (tab[p] < tab[l])
          {
            tablev = tab[p];
            tab[p] = tab[l];
            tab[l] = tablev;
          }
          if (tab[l] < tab[r])
          {
            if (tab[p] > tab[r])
            {
              tablev = tab[r];
              tab[r] = tab[l];
              tab[l] = tablev;
            }
            else
            {
              tablev = tab[p];
              tab[p] = tab[l];
              tab[l] = tablev;
            }
          }
          else
          {
            tablev = tab[l];
          }

          while (1)
          {
            i++;
            while (tab[i] < tablev)
              i++;
            j--;
            while (tab[j] > tablev)
              j--;
            if (j <= i) break;
            tableu = tab[i];
            tab[i] = tab[j];
            tab[j] = tableu;
          }
          tab[l] = tab[j];
          tab[j] = tablev;

          if (r - j > j - l)
          {
            if (j - l > MINI)
            {
              LISTE_L[pstack] = l;
              LISTE_R[pstack] = j - 1;
              pstack++;
              if (pstack >= LSTACK) messageAbort("Stack Overflow");
            }
            l = j + 1;
          }
          else
          {
            if (r - j > MINI)
            {
              LISTE_L[pstack] = j + 1;
              LISTE_R[pstack] = r;
              pstack++;
              if (pstack >= LSTACK) messageAbort("Stack Overflow");
            }
            r = j - 1;
          }
        }
        pstack--;
      }
    }
  }

  /* Final sorting */

  for (j = 1; j < nech; j++)
  {
    if (tab[j - 1] > tab[j])
    {
      tablev = tab[j];
      if (ind) inddev = ind[j];
      i = j;
      while (i > 0)
      {
        if (tab[i - 1] <= tablev) break;
        tab[i] = tab[i - 1];
        if (ind) ind[i] = ind[i - 1];
        i--;
      }
      tab[i] = tablev;
      if (ind) ind[i] = inddev;
    }
  }

  if (safe) tab = (int *) mem_free((char * ) tab);
  return;
}

/****************************************************************************/
/*!
 **  Returns the statistics of an array
 **
 ** \param[in]  nech    Number of samples
 ** \param[in]  tab     Array of values
 ** \param[in]  sel     Array containing the Selection or NULL
 ** \param[in]  wgt     Array containing the Weights or NULL
 **
 ** \param[out]  nval   Number of active values
 ** \param[out]  mini   Minimum value
 ** \param[out]  maxi   Maximum value
 ** \param[out]  delta  Distance between the minimum and the maximum
 ** \param[out]  mean   Mean
 ** \param[out]  stdv   Standard Deviation
 **
 ****************************************************************************/
GEOSLIB_API void ut_statistics(int nech,
                               double *tab,
                               double *sel,
                               double *wgt,
                               int *nval,
                               double *mini,
                               double *maxi,
                               double *delta,
                               double *mean,
                               double *stdv)
{
  int i;
  double num, tmin, tmax, mm, vv, weight;

  /* Initializations */

  tmin = 1.e30;
  tmax = -1.e30;
  num = mm = vv = 0.;
  (*nval) = 0;

  for (i = 0; i < nech; i++)
  {
    if (sel != nullptr && sel[i] == 0.) continue;
    if (FFFF(tab[i])) continue;
    weight = (wgt != nullptr && wgt[i] >= 0) ? wgt[i] :
                                                       1.;
    if (tab[i] < tmin) tmin = tab[i];
    if (tab[i] > tmax) tmax = tab[i];
    (*nval)++;
    num += weight;
    mm += weight * tab[i];
    vv += weight * tab[i] * tab[i];
  }

  /* Returning arguments */

  if (tmax < tmin || (*nval) <= 0)
  {
    *mini = *maxi = *delta = *mean = *stdv = TEST;

  }
  else
  {
    *mini = tmin;
    *maxi = tmax;
    *delta = tmax - tmin;
    mm /= num;
    vv = vv / num - mm * mm;
    if (vv < 0.) vv = 0.;
    *mean = mm;
    *stdv = sqrt(vv);
  }
  return;
}

/****************************************************************************/
/*!
 **  Returns the minimum and maximum of an array
 **
 ** \param[in]  nech    Number of samples
 ** \param[in]  tab     Array of values
 ** \param[in]  sel     Array containing the Selection or NULL
 **
 ** \param[in,out]  nvalid Number of valid samples
 ** \param[in,out]  mini   Minimum value
 ** \param[in,out]  maxi   Maximum value
 **
 ** \remark The calculation starts with the initial values for minimum,
 ** \remark maximum andl nvalid arguments (used in input).
 ** \remark If no default value is defined, set them minimum and maximum to TEST
 ** \remark and set nvalid to 0
 **
 ****************************************************************************/
GEOSLIB_API void ut_stats_mima(int nech,
                               double *tab,
                               double *sel,
                               int *nvalid,
                               double *mini,
                               double *maxi)
{
  double tmin, tmax;
  int i;

  /* Initializations */

  tmin = (FFFF(*mini)) ? 1.e30 :
                         *mini;
  tmax = (FFFF(*maxi)) ? -1.e30 :
                         *maxi;

  for (i = 0; i < nech; i++)
  {
    if (sel != nullptr && sel[i] == 0.) continue;
    if (FFFF(tab[i])) continue;
    if (std::isnan(tab[i])) continue;
    if (std::isinf(tab[i])) continue;
    if (tab[i] < tmin) tmin = tab[i];
    if (tab[i] > tmax) tmax = tab[i];
    (*nvalid)++;
  }

  /* Returning arguments */

  if (tmax < tmin)
  {
    *mini = *maxi = TEST;
  }
  else
  {
    *mini = tmin;
    *maxi = tmax;
  }
  return;
}

/****************************************************************************/
/*!
 **  Print minimum and maximum of an array
 **
 ** \param[in]  title   Title
 ** \param[in]  nech    Number of samples
 ** \param[in]  tab     Array of values
 ** \param[in]  sel     Array containing the Selection or NULL
 **
 ****************************************************************************/
GEOSLIB_API void ut_stats_mima_print(const char *title,
                                     int nech,
                                     double *tab,
                                     double *sel)
{
  int nvalid;
  double mini, maxi;

  // Calculate the min-max statistics

  nvalid = 0;
  mini = maxi = TEST;
  ut_stats_mima(nech, tab, sel, &nvalid, &mini, &maxi);

  // Print the statistics out

  if (nvalid <= 0)
    message("%s: NVal=%6d/%6d - Min=NA - Max=NA\n", title, nvalid, nech);
  else
    message("%s: NVal=%6d/%6d - Min=%lf - Max=%lf\n", title, nvalid, nech, mini,
            maxi);
}

/****************************************************************************/
/*!
 **  Returns the statistics of an array containing the facies
 **
 ** \param[in]  nech    Number of samples
 ** \param[in]  tab     Array of values
 ** \param[in]  sel     Array containing the Selection or NULL
 **
 ** \param[out]  nval   Number of active values
 ** \param[out]  mini   Minimum value
 ** \param[out]  maxi   Maximum value
 **
 ****************************************************************************/
GEOSLIB_API void ut_facies_statistics(int nech,
                                      double *tab,
                                      double *sel,
                                      int *nval,
                                      int *mini,
                                      int *maxi)
{
  int i, number, facies, facmin, facmax;

  /* Initializations */

  facmin = 9999999;
  facmax = 0;
  number = 0;

  for (i = 0; i < nech; i++)
  {
    if (sel != nullptr && sel[i] == 0.) continue;
    if (FFFF(tab[i])) continue;
    facies = (int) tab[i];
    if (facies < 0) continue;
    if (facies < facmin) facmin = facies;
    if (facies > facmax) facmax = facies;
    number++;
  }

  /* Returning arguments */

  if (facmax < facmin || number <= 0)
  {
    *mini = *maxi = ITEST;
  }
  else
  {
    *mini = facmin;
    *maxi = facmax;
    *nval = number;
  }
  return;
}

/****************************************************************************/
/*!
 **  Classify the samples into integer sieves
 **
 ** \param[in]  nech   Number of samples
 ** \param[in]  tab    Array of values
 ** \param[in]  sel    Array containing the Selection or NULL
 ** \param[in]  nclass Number of sieve classes
 ** \param[in]  start  Starting sieve value
 ** \param[in]  pas    Width of the sieve
 **
 ** \param[out]  nmask  Number of masked values
 ** \param[out]  ntest  Number of undefined values
 ** \param[out]  nout   Number of values outside the classes
 ** \param[out]  classe Array for number of samples per sieve
 **
 *****************************************************************************/
GEOSLIB_API void ut_classify(int nech,
                             double *tab,
                             double *sel,
                             int nclass,
                             double start,
                             double pas,
                             int *nmask,
                             int *ntest,
                             int *nout,
                             int *classe)
{
  int i, icl, rank;

  double value;

  /* Initializations */

  for (icl = 0; icl < nclass; icl++)
    classe[icl] = 0;
  (*ntest) = (*nmask) = (*nout) = 0;

  /* Loop on the active information */

  for (i = 0; i < nech; i++)
  {
    if (sel != nullptr && sel[i] == 0.)
    {
      (*nmask)++;
      continue;
    }
    value = tab[i];
    if (FFFF(value))
    {
      (*ntest)++;
      continue;
    }
    rank = (int) ((value - start) / pas);
    if (rank < 0 || rank >= nclass)
    {
      (*nout)++;
      continue;
    }
    classe[rank]++;
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculates the trigonometric features
 **
 ** \param[in]  angle input angle (in degrees)
 **
 ** \param[out]  cosa  cosine function
 ** \param[out]  sina  sine function
 **
 *****************************************************************************/
GEOSLIB_API void ut_rotation_sincos(double angle, double *cosa, double *sina)
{
  double value;

  if (angle == 0.)
  {
    *cosa = 1.;
    *sina = 0.;
    return;
  }
  else if (angle == 90.)
  {
    *cosa = 0.;
    *sina = 1.;
    return;
  }
  else if (angle == 180.)
  {
    *cosa = -1.;
    *sina = 0.;
    return;
  }
  else if (angle == 270.)
  {
    *cosa = 0.;
    *sina = -1.;
    return;
  }
  else
  {
    value = ut_deg2rad(angle);
    *cosa = cos(value);
    *sina = sin(value);
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculates the 2-D rotation matrix
 **
 ** \param[in]  angle Rotation angle (in degrees)
 **
 ** \param[out]  rot   Rotation matrix (Dimension = 4)
 **
 *****************************************************************************/
GEOSLIB_API void ut_rotation_matrix_2D(double angle, double *rot)
{
  double ca, sa;

  ut_rotation_sincos(angle, &ca, &sa);

  /* Define the 2-D rotation matrix */

  rot[0] = ca;
  rot[1] = sa;
  rot[2] = -sa;
  rot[3] = ca;

  return;
}

/*****************************************************************************/
/*!
 **  Calculates the 3-D rotation matrix
 **
 ** \param[in]  alpha angle (in degrees) / oz
 ** \param[in]  beta  angle (in degrees) / oy'
 ** \param[in]  gamma angle (in degrees) / ox''
 **
 ** \param[out] rot   direct rotation matrix (Dimension = 9)
 **
 *****************************************************************************/
GEOSLIB_API void ut_rotation_matrix_3D(double alpha,
                                       double beta,
                                       double gamma,
                                       double *rot)
{
  double ca[3], sa[3];

  /* Initializations */

  ut_rotation_sincos(alpha, &ca[0], &sa[0]);
  ut_rotation_sincos(beta,  &ca[1], &sa[1]);
  ut_rotation_sincos(gamma, &ca[2], &sa[2]);

  /* Define the 3-D rotation matrix */

  rot[0] =  ca[0] * ca[1];
  rot[3] = -sa[0] * ca[2] + ca[0] * sa[1] * sa[2];
  rot[6] =  sa[0] * sa[2] + ca[0] * sa[1] * ca[2];
  rot[1] =  sa[0] * ca[1];
  rot[4] =  ca[0] * ca[2] + sa[0] * sa[1] * sa[2];
  rot[7] = -ca[0] * sa[2] + sa[0] * sa[1] * ca[2];
  rot[2] = -sa[1];
  rot[5] =  ca[1] * sa[2];
  rot[8] =  ca[1] * ca[2];

  return;
}

/*****************************************************************************/
/*!
 **  Calculates the rotation matrix
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  angles Array of angles
 **
 ** \param[out] rot   direct rotation matrix (Dimension = 9)
 **
 *****************************************************************************/
GEOSLIB_API void ut_rotation_matrix(int ndim, const double *angles, double *rot)
{
  if (ndim == 2)
    ut_rotation_matrix_2D(angles[0], rot);
  else if (ndim == 3)
    ut_rotation_matrix_3D(angles[0], angles[1], angles[2], rot);
  else
    ut_rotation_init(ndim, rot);
}

/*****************************************************************************/
/*!
 **  Calculates the rotation matrix.
 **  Returns the rotation matrix as a VectorDouble
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  angles Array of angles
 **
 *****************************************************************************/
GEOSLIB_API VectorDouble ut_rotation_matrix_VD(int ndim,
                                               const VectorDouble& angles)
{
  VectorDouble rot;

  rot.resize(ndim * ndim);
  if (ndim == 2)
    ut_rotation_matrix_2D(angles[0], rot.data());
  else if (ndim == 3)
    ut_rotation_matrix_3D(angles[0], angles[1], angles[2], rot.data());
  else
    ut_rotation_init(ndim, rot.data());

  return rot;
}

/*****************************************************************************/
/*!
 **  Copy a rotation matrix
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  rotin  Input rotation matrix
 **
 ** \param[out] rotout Output rotation matrix (already allocated)
 **
 *****************************************************************************/
GEOSLIB_API void ut_rotation_copy(int ndim, const double *rotin, double *rotout)
{
  int i;

  for (i = 0; i < ndim * ndim; i++)
    rotout[i] = rotin[i];
}

/****************************************************************************/
/*!
 **  Merge the extensions of the boxes (parallel to main axes)
 **
 ** \return  Returns the field extension
 **
 ** \param[in]  ndim     Space dimension
 ** \param[in]  mini1    Input array containing the minimum along each axis
 ** \param[in]  maxi1    Input array containing the maximum along each axis
 ** \param[in]  mini2    Output array containing the minimum along each axis
 ** \param[in]  maxi2    Output array containing the maximum along each axis
 **
 ** \remark  The extension is calculated by merging the extensions of the
 ** \remark  input and output Db structures
 **
 *****************************************************************************/
GEOSLIB_API double ut_merge_extension(int ndim,
                                      double *mini1,
                                      double *maxi1,
                                      double *mini2,
                                      double *maxi2)
{
  double delta, field, mini, maxi;
  int idim;

  /* Loop on the dimensions */

  field = 0.;
  for (idim = 0; idim < ndim; idim++)
  {
    mini = 1.e30;
    if (mini1 != nullptr) mini = MIN(mini, mini1[idim]);
    if (mini2 != nullptr) mini = MIN(mini, mini2[idim]);

    maxi = -1.e30;
    if (maxi1 != nullptr) maxi = MAX(maxi, maxi1[idim]);
    if (maxi2 != nullptr) maxi = MAX(maxi, maxi2[idim]);

    if (mini > maxi) return (TEST);
    delta = maxi - mini;
    field += delta * delta;
  }
  return (sqrt(field));
}

/****************************************************************************/
/*!
 **  Reset the DEBUG status to Idle value
 **
 *****************************************************************************/
GEOSLIB_API void debug_reset(void)

{
  int i;

  for (i = 0; i < DBG_NUMBER; i++)
    DBG[i].status = 0;
  DBGENV.reference = 0;
  DBGENV.index = 0;

  return;
}

/****************************************************************************/
/*!
 **  Print the status of the debug options
 **
 *****************************************************************************/
GEOSLIB_API void debug_print(void)

{
  int i;
  const char *STATUS[] = { "OFF", "ON" };

  mestitle(1, "Parameters for DEBUG option");
  for (i = 0; i < DBG_NUMBER; i++)
    message(". %-30s [%-10s] = %s\n", DBG[i].comment, DBG[i].keyword,
            STATUS[DBG[i].status]);
  if (DBGENV.reference > 0)
    message("Index of the reference target under DEBUG = %d\n",
            DBGENV.reference);
  message(
      "Use 'debug.define' or 'debug.reference' to modify previous values\n");
  return;
}

/****************************************************************************/
/*!
 **  Set the rank of the current target for checking the DEBUG status
 **
 ** \param[in]  rank    rank of the DEBUG target or 0 for undefine
 **
 *****************************************************************************/
GEOSLIB_API void debug_index(int rank)
{

  DBGENV.index = rank;

  return;
}

/****************************************************************************/
/*!
 **  Check if a DEBUG reference index has been defined
 **
 ** \return 0 if DEBUG reference is not defined; rank of this index otherwise
 **
 *****************************************************************************/
GEOSLIB_API int is_debug_reference_defined(void)

{
  return (DBGENV.reference);
}

/****************************************************************************/
/*!
 **  Set the rank of the DEBUG target for an environment variable
 **
 ** \param[in]  rank    rank of the DEBUG target or 0 for undefine
 **
 *****************************************************************************/
GEOSLIB_API void debug_reference(int rank)

{
  DBGENV.reference = rank;

  return;
}

/****************************************************************************/
/*!
 **  Force the action according to the Target Debugging option
 **
 *****************************************************************************/
GEOSLIB_API int debug_force(void)

{
  if (DBGENV.reference <= 0) return (0);
  if (DBGENV.index != DBGENV.reference) return (0);
  return (1);
}

/****************************************************************************/
/*!
 **  Set the DEBUG status for an environment variable
 **
 ** \param[in]  name    name of the environment where DEBUG status is set
 ** \param[in]  status  value of the DEBUG status
 **
 *****************************************************************************/
GEOSLIB_API void debug_define(const char *name, int status)
{
  int i, found;

  /* Look for the "all" keyword */

  if (!strcmp(name, "all"))
  {
    /* Set all the flags to 1 (except the one for Memory) */
    for (i = 0; i < DBG_NUMBER; i++)
      DBG[i].status = status;
    return;
  }

  /* Look for the environment variable */

  for (i = 0, found = -1; i < DBG_NUMBER; i++)
    if (!strcmp(name, DBG[i].keyword)) found = i;

  if (found < 0)
  {
    message("The keywords for DEBUG definition are:\n");
    for (i = 0; i < DBG_NUMBER; i++)
      message("%-10s : %-30s\n", DBG[i].keyword, DBG[i].comment);
    message("The Keyword '%s' is unknown\n", name);
  }
  else
    DBG[found].status = status;

  return;
}

/****************************************************************************/
/*!
 **  Returns the DEBUG status for an environment variable
 **
 ** \return  Debug status
 **
 ** \param[in]  name  name of the environment where DEBUG status is set
 **
 *****************************************************************************/
GEOSLIB_API int debug_query(const char *name)

{
  int i;

  /* The current index coincides with the reference */

  if (debug_force()) return (1);

  /* Look for the environment variable */

  for (i = 0; i < DBG_NUMBER; i++)
    if (!strcmp(name, DBG[i].keyword)) return (DBG[i].status);

  message("The keywords for DEBUG definition are:\n");
  for (i = 0; i < DBG_NUMBER; i++)
    message("%-10s : %-30s\n", DBG[i].keyword, DBG[i].comment);
  message("The Keyword '%s' is unknown\n", name);
  return (0);
}

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
GEOSLIB_API void projec_toggle(int mode)
{
  int projec_actif;

  /* Process the toggling */

  projec_actif = PROJEC.actif;
  if (mode == 1)
    projec_actif = 1;
  else if (mode == 0)
    projec_actif = 0;
  else if (mode == -1) projec_actif = 1 - projec_actif;

  /* Check that no Variety is defined */

  if (VARIETY.flag_sphere && projec_actif)
  {
    messerr("Error when toggling a Projection ON");
    messerr(
        "Definition of a Projection is incompatible with working on a Variety");
    messerr(
        "Cancel the (spherical) Variety first and define the Projection again");
  }
  else
    PROJEC.actif = projec_actif;

  return;
}

/****************************************************************************/
/*!
 **  Toggle the status of the Variety flag
 **
 ** \param[in]  mode Toggle of the Variety flag
 ** \li               0   : Switch the flag OFF
 ** \li               1   : Switch the flag ON
 ** \li              else : Toggle the flag
 **
 *****************************************************************************/
GEOSLIB_API void variety_toggle(int mode)
{
  int variety_actif;

  /* Process the toggling */

  variety_actif = VARIETY.flag_sphere;
  if (mode == 1)
    variety_actif = 1;
  else if (mode == 0)
    variety_actif = 0;
  else if (mode == -1) variety_actif = 1 - variety_actif;

  /* Check that no Variety is defined */

  if (PROJEC.actif && variety_actif)
  {
    messerr("Error when toggling a Spherical Variety ON");
    messerr(
        "Definition of a Variety is incompatible with working on a Projection");
    messerr("Cancel the Projection first and define the Variety again");
  }
  else
    VARIETY.flag_sphere = variety_actif;

  return;
}

/****************************************************************************/
/*!
 **  Returns the projection characteristics
 **
 ** \param[out]  actif activity flag
 **
 *****************************************************************************/
GEOSLIB_API void projec_query(int *actif)

{
  *actif = PROJEC.actif;

  return;
}

/****************************************************************************/
/*!
 **  Print the characteristics of the projection
 **
 *****************************************************************************/
GEOSLIB_API void projec_print(void)

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
 **  Define the Variety characteristics
 **
 ** \param[in]  flag_sphere  1 if the Spherical Variety must be used
 ** \param[in]  radius       Radius of the Sphere
 **
 *****************************************************************************/
GEOSLIB_API void variety_define(int flag_sphere, double radius)
{
  int projec_actif;

  /* Check that no Projection is defined */

  projec_query(&projec_actif);
  if (IFFFF(flag_sphere)) flag_sphere = 0;
  if (flag_sphere && projec_actif)
  {
    messerr("Error when defining a Variety");
    messerr("The definition of a Variety is incompatible with Projections");
    messerr("Cancel the Projection first and define the Variety again");
    return;
  }

  VARIETY.flag_sphere = flag_sphere;
  VARIETY.radius = radius;
  return;
}

/****************************************************************************/
/*!
 **  Returns the Variety presence
 **
 ** \param[out]  flag_sphere 1 if the Spherical coordinates must be used
 **
 *****************************************************************************/
GEOSLIB_API void variety_query(int *flag_sphere)

{
  *flag_sphere = VARIETY.flag_sphere;

  return;
}

/****************************************************************************/
/*!
 **  Returns the Variety characteristics
 **
 ** \param[out]  radius  Radius of the Sphere for the Spherical System
 **
 *****************************************************************************/
GEOSLIB_API void variety_get_characteristics(double *radius)

{
  *radius = VARIETY.radius;

  return;
}

/****************************************************************************/
/*!
 **  Print the characteristics of the Variety
 **
 *****************************************************************************/
GEOSLIB_API void variety_print(void)

{
  if (!VARIETY.flag_sphere) return;
  mestitle(1, "Parameters for Variety Definition");
  message("The Spherical Variety is defined\n");
  message("- Radius of the Sphere = %lf\n", VARIETY.radius);
  return;
}


/****************************************************************************/
/*!
 **  Ask for the characteristics of a matrix
 **
 ** \param[in]  title     Title of the matrix
 ** \param[in]  flag_sym  1 if the matrix MUST be symmetrical
 ** \param[in]  flag_def  1 for default values; 0 otherwise
 ** \param[in]  nx        Number of columns
 ** \param[in]  ny        Number of rows
 ** \param[in]  valmin    Minimum authorized value or TEST
 ** \param[in]  valmax    Maximum authorized value or TEST
 ** \param[in,out] tab    Input/Output matrix (if flag_def=1)
 **
 *****************************************************************************/
GEOSLIB_API void get_matrix(const char *title,
                            int flag_sym,
                            int flag_def,
                            int nx,
                            int ny,
                            double valmin,
                            double valmax,
                            double *tab)
{
  int ix, iy;

  /* Loop on the matrix elements */

  for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++)
    {
      if (flag_sym && iy > ix) continue;
      (void) gslSPrintf(QUESTION, "%s(", title);
      if (nx > 1) (void) gslSPrintf(&QUESTION[strlen(QUESTION)], "%d", ix + 1);
      if (nx > 1 && ny > 1) (void) gslStrcat(QUESTION, ",");
      if (ny > 1) (void) gslSPrintf(&QUESTION[strlen(QUESTION)], "%d", iy + 1);
      (void) gslStrcat(QUESTION, ")");
      TAB(ix,iy) = _lire_double(QUESTION, flag_def, (flag_def) ? TAB(ix, iy) :
                                                                 TEST,
                                valmin, valmax);
      if (flag_sym) TAB(iy,ix) = TAB(ix, iy);
    }

  return;
}

/****************************************************************************/
/*!
 **  Ask for the characteristics of a rotation
 **
 ** \param[in]  title     Title of the rotation
 ** \param[in]  flag_def  1 for default values; 0 otherwise
 ** \param[in]  ndim      Dimension of the space
 ** \param[in,out] rot    Input/Output rotation matrix (if flag_def=1)
 **
 *****************************************************************************/
GEOSLIB_API void get_rotation(const char *title,
                              int flag_def,
                              int ndim,
                              double *rot)
{
  double dir[2], angles2D[2], angles3D[3], alpha, beta, gamma;
  int mode;

  /* Print the title */

  mode = 1;
  message("%s : \n", title);

  /* Choose the rotation mode according to the space dimension */

  switch (ndim)
  {
    case 1:
      messerr("The rotation should not be defined in 1-D space");
      break;

    case 2:
      message("Rotation Definition mode :\n");
      message("1 - The rotation angle(s) in degrees\n");
      message("2 - The rotation matrix\n");
      message("3 - The vector defining the X-axis after rotation\n");
      mode = _lire_int("Choose the Definition mode", 1, 1, 1, 3);
      break;

    case 3:
      message("Rotation Definition mode :\n");
      message("1 - The rotation angle(s) in degrees\n");
      message("2 - The rotation matrix\n");
      mode = _lire_int("Choose the Definition mode", 1, 1, 1, 3);
      break;
  }

  /* Define the rotation */

  switch (mode)
  {
    case 1: /* Rotation angles */
      if (ndim == 2)
      {
        (void) ut_angles_from_rotation_matrix(rot, ndim, angles2D);
        alpha = _lire_double("Rotation angle from East counterclockwise", 1,
                             angles2D[0], TEST, TEST);
        ut_rotation_matrix_2D(alpha, rot);
      }
      else
      {
        (void) ut_angles_from_rotation_matrix(rot, ndim, angles3D);
        alpha = _lire_double("Rotation angle around Oz  ", 1, angles3D[0], TEST,
                             TEST);
        beta = _lire_double("Rotation angle around Oy' ", 1, angles3D[1], TEST,
                            TEST);
        gamma = _lire_double("Rotation angle around Ox''", 1, angles3D[2], TEST,
                             TEST);
        ut_rotation_matrix_3D(alpha, beta, gamma, rot);
      }
      break;

    case 2: /* Rotation matrix */
      get_matrix("Rotation Matrix", 0, flag_def, ndim, ndim, -1., 1., rot);
      break;

    case 3: /* The rotation vector (only for 2-D) */
      dir[0] = 1.;
      dir[1] = 0.;
      get_matrix("Rotation Vector", 0, flag_def, 2, 1, -1., 1., dir);
      ROT(0,0) = dir[0];
      ROT(1,0) = dir[1];
      ROT(0,1) = -dir[1];
      ROT(1,1) = dir[0];
      break;
  }

  return;
}

/****************************************************************************/
/*!
 **  Calculates the rotation angles from the rotation matrix
 **
 ** \param[in]  rot   Rotation matrix (Dimension = 9)
 ** \param[in]  ndim  Space dimension
 **
 ** \param[out]  angles Rotation angles (Dimension = ndim)
 **
 *****************************************************************************/
GEOSLIB_API int ut_angles_from_rotation_matrix(const double *rot,
                                               int ndim,
                                               double *angles)
{
  double s0, c0, s1, c1, s2, c2;
  int i, nval;

  /* Initializations */

  for (i = 0; i < ndim; i++)
    angles[i] = 0.;
  if (rot == nullptr) return (0);

  /* Dispatch */

  if (ndim == 1)
  {
    nval = 1;
  }
  else if (ndim == 2)
  {
    angles[0] = atan2(rot[1], rot[0]);
    nval = 1;
  }
  else if (ndim == 3)
  {
    nval = 3;
    s1 = -rot[2];
    c1 = sqrt(rot[0] * rot[0] + rot[1] * rot[1]);
    if (ABS(c1) < EPSILON10)
    {
      if (s1 > 0.)
      {
        angles[0] = 0.;
        angles[1] = GV_PI / 2.;
        angles[2] = atan2(rot[3], rot[6]);
      }
      else
      {
        angles[0] = 0.;
        angles[1] = -GV_PI / 2.;
        angles[2] = atan2(-rot[3], -rot[6]);
      }
    }
    else
    {
      c0 = rot[0] / c1;
      s0 = rot[1] / c1;
      s2 = rot[5] / c1;
      c2 = rot[8] / c1;
      angles[0] = atan2(s0, c0);
      angles[1] = atan2(s1, c1);
      angles[2] = atan2(s2, c2);
    }
  }
  else
  {
    nval = 0;
  }

  /* Convert into degrees */

  for (i = 0; i < nval; i++)
    angles[i] = ut_rad2deg(angles[i]);

  return (nval);
}

/****************************************************************************/
/*!
 **  Calculates the rotation angle from the direction coefficient
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  ndir   Number of directions
 ** \param[in]  codir  Direction vector (Dimension = ndim * ndir)
 **
 ** \param[out]  angles Rotation angles (Dimension = ndim * ndir)
 **
 *****************************************************************************/
GEOSLIB_API void ut_angles_from_codir(int ndim,
                                      int ndir,
                                      const VectorDouble& codir,
                                      VectorDouble& angles)
{
  double norme;
  int i, nval;

  /* Initializations */

  for (i = 0; i < ndim * ndir; i++) angles[i] = 0.;

  /* Dispatch */

  if (ndim == 1)
  {
    return;
  }
  else if (ndim == 2)
  {
    angles[0] = atan2(codir[1], codir[0]);
    angles[1] = 0.;
    nval = 1;
  }
  else if (ndim == 3)
  {
    norme = codir[0] * codir[0] + codir[1] * codir[1];
    if (norme > 0.)
    {
      norme = sqrt(norme);
      angles[0] = atan2(codir[1] / norme, codir[0] / norme);
      angles[1] = atan2(codir[2], norme);
    }
    nval = 2;
  }
  else
  {
    nval = 0;
  }

  /* Convert into degrees */

  for (i = 0; i < nval; i++)
    angles[i] = ut_rad2deg(angles[i]);

  return;
}

/****************************************************************************/
/*!
 **  Starting from a rotation matrix, check it is different from the Identity
 **
 ** \return  1 if a rotation is defined; 0 otherwise
 **
 ** \param[in]  rot      Rotation matrix
 ** \param[in]  ndim     Space dimension
 **
 *****************************************************************************/
GEOSLIB_API int ut_rotation_check(double *rot, int ndim)
{
  int i, j;

  for (i = 0; i < ndim; i++)
    for (j = 0; j < ndim; j++)
    {
      if (i == j)
      {
        if (ABS(ROT(i,j) - 1.) > EPSILON10) return (1);
      }
      else
      {
        if (ABS(ROT(i,j)) > EPSILON10) return (1);
      }
    }
  return (0);
}

/****************************************************************************/
/*!
 **  Golden Search algorithm
 **
 ** \return  The estimated value
 **
 ** \param[in]  func_evaluate  Evaluating function
 ** \param[in]  user_data      User Data
 ** \param[in]  tolstop        Tolerance parameter
 ** \param[in]  a0             Initial value for lower bound of interval
 ** \param[in]  c0             Initial value for upper bound of interval
 **
 ** \param[out] test_loc       Final value of the evaluating function
 ** \param[out] niter          Number of iterations
 **
 *****************************************************************************/
GEOSLIB_API double golden_search(double (*func_evaluate)(double test,
                                                         void *user_data),
                                 void *user_data,
                                 double tolstop,
                                 double a0,
                                 double c0,
                                 double *test_loc,
                                 double *niter)
{
  double phi, resphi, b, x, fb, fx, result, a, c;
  int flag_test;

  /* Initializations */

  phi = (1. + sqrt(5.)) / 2.;
  resphi = 2. - phi;
  a = a0;
  c = c0;

  /* Initial values for the golden search */

  b = (a + c) / 2;
  fb = func_evaluate(b, user_data);

  (*niter) = 1;
  while ((c - a) > tolstop)
  {
    flag_test = (c - b) > (b - a);
    if (flag_test)
      x = b + resphi * (c - b);
    else
      x = b - resphi * (b - a);
    fx = func_evaluate(x, user_data);
    (*niter) = (*niter) + 1.;

    if (fx < fb)
    {
      if (flag_test)
      {
        a = b;
        b = x;
        fb = fx;
      }
      else
      {
        c = b;
        b = x;
        fb = fx;
      }
    }
    else
    {
      if (flag_test)
        c = x;
      else
        a = x;
    }
  }

  /* Returning value */

  result = (a + c) / 2.;
  if (test_loc != NULL) *test_loc = fb;

  return (result);
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
static Keypair *st_get_keypair_address(const char *keyword)

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
    KEYPAIR_TABS = (Keypair *)
        realloc((char *) KEYPAIR_TABS, sizeof(Keypair) * KEYPAIR_NTAB);
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
        free((char *) keypair->values);
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
        messageAbort("Keypair append cannot change origin or number of columns");
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

    keypair->values = (double *) malloc(sizeof(double) * new_size);
  }
  else
  {

    // The old_dimensions are non zero, reallocate the contents

    keypair->values = (double *) realloc((char *) keypair->values,
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
    icopy = (int *) values;
    for (int i = 0; i < size; i++)
      keypair->values[i + start] = icopy[i];
  }
  else
  {
    rcopy = (double *) values;
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
GEOSLIB_API void set_keypair(const char *keyword,
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

  st_keypair_copy(keypair, 2, 0, (void *) values);

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
GEOSLIB_API void app_keypair(const char *keyword,
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

  st_keypair_copy(keypair, 2, start, (void *) values);
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
GEOSLIB_API void set_keypair_int(const char *keyword,
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

  st_keypair_copy(keypair, 1, 0, (void *) values);
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
GEOSLIB_API void app_keypair_int(const char *keyword,
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

  st_keypair_copy(keypair, 1, start, (void *) values);
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
  free((char *) keypair->values);
  keypair->values = nullptr;

  /* Shift all subsequent keypairs */

  for (int i = indice + 1; i < KEYPAIR_NTAB; i++)
    KEYPAIR_TABS[i - 1] = KEYPAIR_TABS[i];

  KEYPAIR_NTAB--;
  KEYPAIR_TABS = (Keypair *) realloc((char *) KEYPAIR_TABS,
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
GEOSLIB_API void del_keypair(const char *keyword, int flag_exact)
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
GEOSLIB_API double get_keypone(const char *keyword, double valdef)
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
    rtab = (double *) keypair->values;
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
GEOSLIB_API int get_keypair(const char *keyword,
                            int *nrow,
                            int *ncol,
                            double **values)
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

  valloc = (double *) malloc(sizeof(double) * size);
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
GEOSLIB_API int get_keypair_int(const char *keyword,
                                int *nrow,
                                int *ncol,
                                int **values)
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

  valloc = (int *) malloc(sizeof(int) * size);
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
GEOSLIB_API void print_keypair(int flag_short)

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
 **  Initialize a rotation matrix
 **
 ** \param[in]  ndim      Space dimension
 **
 ** \param[out] rot       Rotation matrix
 **
 *****************************************************************************/
GEOSLIB_API void ut_rotation_init(int ndim, double *rot)
{
  int i, j, ecr;

  for (i = ecr = 0; i < ndim; i++)
    for (j = 0; j < ndim; j++, ecr++)
      rot[ecr] = (i == j);
}

/*****************************************************************************/
/*!
 **  Calculate the median from a table of values
 **
 ** \returns The median value
 **
 ** \param[in]  tab       Array of values
 ** \param[in]  ntab      Number of samples
 **
 *****************************************************************************/
GEOSLIB_API double ut_median(double *tab, int ntab)
{
  int i, j, k, nr, nl, even, lo, hi, loop, mid;
  double result, xlo, xhi, temp, xmin, xmax;

  nr = ntab / 2;
  nl = nr - 1;
  even = 0;
  /* hi & lo are position limits encompassing the median. */
  lo = 0;
  hi = ntab - 1;

  if (ntab == 2 * nr) even = 1;
  if (ntab < 3)
  {
    if (ntab < 1) return 0.;
    if (ntab == 1) return tab[0];
    return 0.5 * (tab[0] + tab[1]);
  }

  /* Find median of 1st, middle & last values. */
  do
  {
    mid = (lo + hi) / 2;
    result = tab[mid];
    xlo = tab[lo];
    xhi = tab[hi];
    if (xhi < xlo)
    {
      temp = xlo;
      xlo = xhi;
      xhi = temp;
    }
    if (result > xhi)
      result = xhi;
    else if (result < xlo) result = xlo;
    /* The basic quicksort algorithm to move all values <= the sort key (XMED)
     * to the left-hand end, and all higher values to the other end.
     */
    i = lo;
    j = hi;
    do
    {
      while (tab[i] < result)
        i++;
      while (tab[j] > result)
        j--;
      loop = 0;
      if (i < j)
      {
        temp = tab[i];
        tab[i] = tab[j];
        tab[j] = temp;
        i++;
        j--;
        if (i <= j) loop = 1;
      }
    }
    while (loop); /* Decide which half the median is in. */

    if (even)
    {
      if (j == nl && i == nr)
      /* Special case, n even, j = n/2 & i = j + 1, so the median is
       * between the two halves of the series.   Find max. of the first
       * half & min. of the second half, then average.
       */
      {
        xmax = tab[0];
        xmin = tab[ntab - 1];
        for (k = lo; k <= j; k++)
          xmax = MAX(xmax, tab[k]);
        for (k = i; k <= hi; k++)
          xmin = MIN(xmin, tab[k]);
        return 0.5 * (xmin + xmax);
      }
      if (j < nl) lo = i;
      if (i > nr) hi = j;
      if (i == j)
      {
        if (i == nl) lo = nl;
        if (j == nr) hi = nr;
      }
    }
    else
    {
      if (j < nr) lo = i;
      if (i > nr) hi = j;
      if (i == j && i == nr) return result;
    }
  }
  while (lo < hi - 1);

  if (even) return (0.5 * (tab[nl] + tab[nr]));
  if (tab[lo] > tab[hi])
  {
    temp = tab[lo];
    tab[lo] = tab[hi];
    tab[hi] = temp;
  }
  return tab[nr];
}

/****************************************************************************/
/*!
 **  Normalize a vector
 **
 ** \param[in]  ntab   Vector dimension
 ** \param[in,out]  tab    Vector to be normalized
 **
 *****************************************************************************/
GEOSLIB_API void ut_normalize(int ntab, double *tab)
{
  int i;
  double norme;

  norme = 0.;
  for (i = 0; i < ntab; i++)
    norme += tab[i] * tab[i];
  norme = sqrt(norme);

  if (norme <= 0.) return;
  for (i = 0; i < ntab; i++)
    tab[i] /= norme;
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
GEOSLIB_API int solve_P2(double a, double b, double c, double *x)
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
GEOSLIB_API int solve_P3(double a, double b, double c, double d, double *x)
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
      u = (u > 0.) ? pow(u, 1. / 3.) :
                     -pow(-u, 1. / 3.);
      v = (-q - s1) / 2.;
      v = (v > 0.) ? pow(v, 1. / 3.) :
                     -pow(-v, 1. / 3.);
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
        x[k] = ecart
            + 2. * sqrt(-p / 3.) * cos((acos(s1) + 2. * k * GV_PI) / 3.);
      return (3);
    }
  }
}

/****************************************************************************/
/*!
 **  Manage the PL_Dist structure
 **
 ** \return Pointer to the PL_Dist structure
 **
 ** \param[in]  mode       Management operation
 ** \li                    1 : Allocation
 ** \li                   -1 : Deallocation
 ** \param[in]  pldist_loc Input PL_Dist structure (used for mode=-1)
 ** \param[in]  ndim       Space dimension
 **
 ** \remarks The PL_Dist structure that has been allocated (mode=1),
 ** \remarks must be freed using the same function with mode=-1
 **
 *****************************************************************************/
GEOSLIB_API PL_Dist *pldist_manage(int mode,
                                   PL_Dist *pldist_loc,
                                   int ndim,
                                   int /*nvert*/)
{
  PL_Dist *pldist;
  int idim;

  /* Dispatch */

  if (mode > 0)
  {
    pldist = (PL_Dist *) mem_alloc(sizeof(PL_Dist), 1);
    pldist->ndim = ndim;
    pldist->rank = -1;
    pldist->dist = TEST;
    pldist->coor = (double *) mem_alloc(sizeof(double) * ndim, 1);
    for (idim = 0; idim < ndim; idim++)
      pldist->coor[idim] = TEST;
  }
  else
  {
    pldist = pldist_loc;
    if (pldist == (PL_Dist *) NULL) return (pldist);
    pldist->coor = (double *) mem_free((char * ) pldist->coor);
    pldist = (PL_Dist *) mem_free((char * ) pldist);
  }

  return (pldist);
}

/****************************************************************************/
/*!
 **  Find the shortest distance between the point (x0,y0) and the segment
 **  with the two end points (x1,y1) and (x2,y2)
 **
 ** \return Minimum algebraic distance (positive or negative)
 **
 ** \param[in]  x0,y0   Coordinates of the target point
 ** \param[in]  x1,y1   Coordinate of the first end-point of the segment
 ** \param[in]  x2,y2   Coordinate of the second end-point of the segment
 **
 ** \param[out] xd,yd   Coordinates of the closest point
 ** \param[out] nint    =1 if the projection belongs to the segment
 **                     =0 if it is set to one of the segment vertices
 **
 *****************************************************************************/
GEOSLIB_API double distance_point_to_segment(double x0,
                                             double y0,
                                             double x1,
                                             double y1,
                                             double x2,
                                             double y2,
                                             double *xd,
                                             double *yd,
                                             int *nint)
{
  double dx, dy, dxp, dyp, ratio, dist, signe;

  dx = x2 - x1;
  dy = y2 - y1;

  ratio = (dx * (x0 - x1) + dy * (y0 - y1)) / (dx * dx + dy * dy);

  if (ratio < 0)
  {
    *xd = x1;
    *yd = y1;
    *nint = 0;
  }
  else if (ratio > 1)
  {
    *xd = x2;
    *yd = y2;
    *nint = 0;
  }
  else
  {
    *xd = x1 + ratio * dx;
    *yd = y1 + ratio * dy;
    *nint = 1;
  }

  dxp = x0 - (*xd);
  dyp = y0 - (*yd);
  dist = sqrt(dxp * dxp + dyp * dyp);
  signe = (dy * (x0 - x1) - dx * (y0 - y1) > 0.) ? 1 :
                                                   -1;

  return (dist * signe);
}

/****************************************************************************/
/*!
 **  Find the shortest distance between the point (x0,y0) and a polyline
 **
 ** \param[in]  x0,y0   Coordinates of the target point
 ** \param[in]  nvert   Number of segments in the polyline
 ** \param[in]  xl      Array of X-coordinates of the polyline
 ** \param[in]  yl      Array of Y-coordinates of the polyline
 **
 ** \param[out] pldist  PL_Dist structure
 **
 ** \remarks  The number of points of the polyline is equal to nvert
 **
 *****************************************************************************/
GEOSLIB_API void distance_point_to_polyline(double x0,
                                            double y0,
                                            int nvert,
                                            const double *xl,
                                            const double *yl,
                                            PL_Dist *pldist)
{
  double xx, yy, dist, dmin;
  int i, nint;

  /* Dispatch */

  dmin = 1.e30;
  for (i = 0; i < nvert - 1; i++)
  {
    dist = distance_point_to_segment(x0, y0, xl[i], yl[i], xl[i + 1], yl[i + 1],
                                     &xx, &yy, &nint);
    if (ABS(dist) > dmin) continue;
    pldist->rank = i;
    pldist->coor[0] = xx;
    pldist->coor[1] = yy;
    pldist->dist = dmin = ABS(dist);
  }
  return;
}

/****************************************************************************/
/*!
 **  Find the shortest distance between two points (x1,y1) and (x2,y2)
 **  which belong to the polyline
 **
 ** \return Minimum distance
 **
 ** \param[in]  pldist1 First PL_Dist structure
 ** \param[in]  pldist2 Second PL_Dist structure
 ** \param[in]  xl      Array of X-coordinates of the polyline
 ** \param[in]  yl      Array of Y-coordinates of the polyline
 **
 ** \remarks  The number of points of the polyline is equal to nvert
 **
 *****************************************************************************/
GEOSLIB_API double distance_along_polyline(PL_Dist *pldist1,
                                           PL_Dist *pldist2,
                                           double *xl,
                                           double *yl)
{
  int i;
  double dist, local1[2], local2[2];
  PL_Dist *pl1, *pl2;

  /* Initializations */

  dist = 0.;
  if (pldist1->rank < pldist2->rank)
  {
    pl1 = pldist1;
    pl2 = pldist2;
  }
  else
  {
    pl1 = pldist2;
    pl2 = pldist1;
  }

  /* If both projected points belong to the same segment */

  if (pl1->rank == pl2->rank)
  {
    dist += ut_distance(2, pl1->coor, pl2->coor);
  }
  else
  {

    /* Distance on the first segment */

    local1[0] = xl[pl1->rank + 1];
    local1[1] = yl[pl1->rank + 1];
    dist += ut_distance(2, pl1->coor, local1);

    /* Distance on the last segment */

    local2[0] = xl[pl2->rank + 1];
    local2[1] = yl[pl2->rank + 1];
    dist += ut_distance(2, pl2->coor, local2);

    for (i = pl1->rank + 1; i < pl2->rank; i++)
    {
      local1[0] = xl[i + 1];
      local1[1] = yl[i + 1];
      local2[0] = xl[i];
      local2[1] = yl[i];
      dist += ut_distance(2, local1, local2);
    }
  }
  return (dist);
}

/****************************************************************************/
/*!
 **  Shift a point along a segment
 **
 ** \param[in]  x1,y1   Coordinates of the first point
 ** \param[in]  x2,y2   Coordinates of the second point
 ** \param[in]  ratio   Shifting ratio
 **
 ** \param[out] x0,y0   Shifted point
 **
 ** \remarks 'ratio' varies between 0 and 1
 ** \remarks When 'ratio' =0, (x0,y0) coincides with (x1,y1)
 ** \remarks When 'ratio'>=1, (x0,y0) coincides with (x2,y2)
 **
 *****************************************************************************/
static void st_shift_point(double x1,
                           double y1,
                           double x2,
                           double y2,
                           double ratio,
                           double *x0,
                           double *y0)
{
  if (ratio <= 0.)
  {
    *x0 = x1;
    *y0 = y1;
  }
  else if (ratio >= 1.)
  {
    *x0 = x2;
    *y0 = y2;
  }
  else
  {
    *x0 = x1 + ratio * (x2 - x1);
    *y0 = y1 + ratio * (y2 - y1);
  }
}

/****************************************************************************/
/*!
 **  Find the shortest distance between two points (x1,y1) and (x2,y2)
 **  passing through a polyline
 **
 ** \return Minimum distance
 **
 ** \param[in]  ap      Coefficient applied to the projected distances
 ** \param[in]  al      Coefficient applied to the distance along line
 ** \param[in]  x1,y1   Coordinates of the first point
 ** \param[in]  x2,y2   Coordinates of the second point
 ** \param[in]  nvert   Number of segments in the polyline
 ** \param[in]  xl      Array of X-coordinates of the polyline
 ** \param[in]  yl      Array of Y-coordinates of the polyline
 **
 ** \remarks  The number of points of the polyline is equal to nvert
 **
 *****************************************************************************/
GEOSLIB_API double distance_points_to_polyline(double ap,
                                               double al,
                                               double x1,
                                               double y1,
                                               double x2,
                                               double y2,
                                               int nvert,
                                               double *xl,
                                               double *yl)
{
  double dist, d1, d2, dh, dv, dloc, dmin, xp1, xp2, yp1, yp2, dist1, dist2;
  PL_Dist *pldist1, *pldist2;

  /* Initialization */

  pldist1 = pldist_manage(1, NULL, 2, nvert);
  pldist2 = pldist_manage(1, NULL, 2, nvert);

  /* Calculate the projection of each end point */

  distance_point_to_polyline(x1, y1, nvert, xl, yl, pldist1);
  distance_point_to_polyline(x2, y2, nvert, xl, yl, pldist2);

  /* Calculate the minimum distance */

  dist = 1.e30;
  dh = dv = 0.;
  d1 = pldist1->dist;
  d2 = pldist2->dist;
  dh = ap * ABS(d1 - d2);

  if (al > 0.)
  {
    dv = distance_along_polyline(pldist1, pldist2, xl, yl);
    d1 = ABS(d1);
    d2 = ABS(d2);
    dmin = MIN(d1, d2);

    xp1 = pldist1->coor[0];
    yp1 = pldist1->coor[1];
    xp2 = pldist2->coor[0];
    yp2 = pldist2->coor[1];
    dist1 = (xp1 - xp2) * (xp1 - xp2) + (yp1 - yp2) * (yp1 - yp2);
    if (ABS(d1) > 0.) st_shift_point(xp1, yp1, x1, y1, dmin / d1, &xp1, &yp1);
    if (ABS(d2) > 0.) st_shift_point(xp2, yp2, x2, y2, dmin / d2, &xp2, &yp2);
    dist2 = (xp1 - xp2) * (xp1 - xp2) + (yp1 - yp2) * (yp1 - yp2);
    dv = (dist1 <= 0.) ? 0. :
                         dv * al * sqrt(dist2 / dist1);
  }
  dloc = sqrt(dh * dh + dv * dv);
  if (dloc < dist) dist = dloc;

  pldist1 = pldist_manage(-1, pldist1, 2, nvert);
  pldist2 = pldist_manage(-1, pldist2, 2, nvert);
  return (dist);
}

/****************************************************************************/
/*!
 **  Convert a string into lowercase characters
 **
 ** \param[in,out]  string  Input/Output string
 **
 *****************************************************************************/
GEOSLIB_API void string_to_lowercase(char *string)

{
  int i, n;

  n = static_cast<int> (strlen(string));
  for (i = 0; i < n; i++)
    if (string[i] >= 'A' && string[i] <= 'Z')
      string[i] = ('a' + string[i] - 'A');
}

/****************************************************************************/
/*!
 **  Convert a string into uppercase characters
 **
 ** \param[in,out]  string  Input/Output string
 **
 *****************************************************************************/
GEOSLIB_API void string_to_uppercase(char *string)

{
  int i, n;

  n = static_cast<int> (strlen(string));
  for (i = 0; i < n; i++)
    if (string[i] >= 'a' && string[i] <= 'z')
      string[i] = ('A' + string[i] - 'a');
}

/****************************************************************************/
/*!
 **  Compare two strings (case dependent or not)
 **
 ** \return 0 if strings are equal; 1 otherwise
 **
 ** \param[in]  flag_case 1 if match must be case independent; 0 otherwise
 ** \param[in]  string1   First input string
 ** \param[in]  string2   Second input string
 **
 *****************************************************************************/
GEOSLIB_API int string_compare(int flag_case,
                               const char *string1,
                               const char *string2)
{
  int flag_diff;

  /* Dispatch */

  if (flag_case)
  {
    flag_diff = strcmp(string1, string2);
  }
  else
  {
    (void) gslStrcpy(INSTR1, string1);
    (void) gslStrcpy(INSTR2, string2);
    string_to_lowercase(INSTR1);
    string_to_lowercase(INSTR2);
    flag_diff = strcmp(INSTR1, INSTR2);
  }
  return (flag_diff);
}

/****************************************************************************/
/*!
 **  Compute combinations(n,k)
 **
 ** \return Return the number of combinations of 'k' objects amongst 'n'
 **
 ** \param[in]  n     Total number of objects (>= 1)
 ** \param[in]  k     Selected number of objects (>= 1)
 **
 *****************************************************************************/
GEOSLIB_API double ut_cnp(int n, int k)
{
  double result, v1, v2;

  result = 0.;
  if (k > n) return (result);

  v1 = v2 = 0.;
  for (int i = 0; i < k; i++)
  {
    v1 += log(n - i);
    v2 += log(i + 1);
  }
  result = exp(v1 - v2);
  return (result);
}

/****************************************************************************/
/*!
 **  Create the matrix containing the Pascal Triangle coefficients
 **
 ** \return A matrix (Dimension: ndim * ndim) containing the coefficients
 ** \return or NULL if core allocation problem has been encountered
 **
 ** \param[in]  ndim   Size of the matrix
 **
 ** \remarks The calling function must free the returned matrix
 **
 *****************************************************************************/
GEOSLIB_API double *ut_pascal(int ndim)
{
  double *m;
#define M(j,i)            (m[(i) * ndim + (j)])

  /* Core allocation */

  m = (double *) mem_alloc(sizeof(double) * ndim * ndim, 0);
  if (m == nullptr) return (m);
  for (int i = 0; i < ndim * ndim; i++)
    m[i] = 0.;

  /* Fill the matrix */

  for (int i = 0; i < ndim; i++)
    for (int j = i; j < ndim; j++)
    {
      if (j == 0 || i == 0)
        M(i,j) = 1.;
      else
        M(i,j) = M(i,j-1) + M(i - 1, j - 1);
    }
  return (m);
#undef M  
}

/****************************************************************************/
/*!
 **  Calculate the geodetic angular distance between two points on the sphere
 **
 ** \return Angular distance
 **
 ** \param[in]  long1  Longitude of the first point (in degrees)
 ** \param[in]  lat1   Latitude of the first point (in degrees)
 ** \param[in]  long2  Longitude of the second point (in degrees)
 ** \param[in]  lat2   Latitude of the second point (in degrees)
 **
 *****************************************************************************/
GEOSLIB_API double ut_geodetic_angular_distance(double long1,
                                                double lat1,
                                                double long2,
                                                double lat2)
{
  double rlon1, rlat1, rlon2, rlat2, dlong, angdst;

  rlon1 = ut_deg2rad(long1);
  rlat1 = ut_deg2rad(lat1);
  rlon2 = ut_deg2rad(long2);
  rlat2 = ut_deg2rad(lat2);
  dlong = rlon2 - rlon1;
  angdst = acos(sin(rlat1) * sin(rlat2) + cos(rlat1) * cos(rlat2) * cos(dlong));
  return (angdst);
}

/****************************************************************************/
/*!
 **  Extract A from a,b,c
 **
 ** \param[in]  cosa   Cosine of first angle
 ** \param[in]  sinb   Sine of second angle
 ** \param[in]  cosb   Cosine of second angle
 ** \param[in]  sinc   Sine of third angle
 ** \param[in]  cosc   Cosine of third angle
 **
 *****************************************************************************/
static double st_convert_geodetic_angle(double /*sina*/,
                                        double cosa,
                                        double sinb,
                                        double cosb,
                                        double sinc,
                                        double cosc)
{
  double prod, cosA;

  prod = sinb * sinc;
  cosA = (prod == 0.) ? 0. : (cosa - cosb * cosc) / prod;
  if (cosA < -1) cosA = -1.;
  if (cosA > +1) cosA = +1.;
  return (acos(cosA));
}

/****************************************************************************/
/*!
 **  Calculate all geodetic angles from a spherical triangle
 **
 ** \param[in]  long1  Longitude of the first point (in degrees)
 ** \param[in]  lat1   Latitude of the first point (in degrees)
 ** \param[in]  long2  Longitude of the second point (in degrees)
 ** \param[in]  lat2   Latitude of the second point (in degrees)
 ** \param[in]  long3  Longitude of the third point (in degrees)
 ** \param[in]  lat3   Latitude of the third point (in degrees)
 **
 ** \param[out] a      Angle (P2,O,P3)
 ** \param[out] b      Angle (P3,O,P1)
 ** \param[out] c      Angle (P1,O,P2)
 ** \param[out] A      Angle (P2,P1,P3)
 ** \param[out] B      Angle (P3,P2,P1)
 ** \param[out] C      Angle (P1,P3,P2)
 **
 *****************************************************************************/
GEOSLIB_API void ut_geodetic_angles(double long1,
                                    double lat1,
                                    double long2,
                                    double lat2,
                                    double long3,
                                    double lat3,
                                    double *a,
                                    double *b,
                                    double *c,
                                    double *A,
                                    double *B,
                                    double *C)
{
  double cosa, cosb, cosc, sina, sinb, sinc;

  *a = ut_geodetic_angular_distance(long2, lat2, long3, lat3);
  *b = ut_geodetic_angular_distance(long1, lat1, long3, lat3);
  *c = ut_geodetic_angular_distance(long1, lat1, long2, lat2);

  cosa = cos(*a);
  cosb = cos(*b);
  cosc = cos(*c);
  sina = sin(*a);
  sinb = sin(*b);
  sinc = sin(*c);

  *A = st_convert_geodetic_angle(sina, cosa, sinb, cosb, sinc, cosc);
  *B = st_convert_geodetic_angle(sinb, cosb, sinc, cosc, sina, cosa);
  *C = st_convert_geodetic_angle(sinc, cosc, sina, cosa, sinb, cosb);

  return;
}

/****************************************************************************/
/*!
 **  Calculate the perimeter of the spherical triangle
 **
 ** \return The Perimeter of the spherical triangle
 **
 ** \param[in]  long1  Longitude of the first point (in degrees)
 ** \param[in]  lat1   Latitude of the first point (in degrees)
 ** \param[in]  long2  Longitude of the second point (in degrees)
 ** \param[in]  lat2   Latitude of the second point (in degrees)
 ** \param[in]  long3  Longitude of the third point (in degrees)
 ** \param[in]  lat3   Latitude of the third point (in degrees)
 **
 *****************************************************************************/
GEOSLIB_API double ut_geodetic_triangle_perimeter(double long1,
                                                  double lat1,
                                                  double long2,
                                                  double lat2,
                                                  double long3,
                                                  double lat3)
{
  double a, b, c, ga, gb, gc, perimeter;

  ut_geodetic_angles(long1, lat1, long2, lat2, long3, lat3, &a, &b, &c, &ga,
                     &gb, &gc);
  perimeter = a + b + c;
  return (perimeter);
}

/****************************************************************************/
/*!
 **  Calculate the surface of the spherical triangle
 **
 ** \return The Surface of the spherical triangle (with unit radius)
 **
 ** \param[in]  long1  Longitude of the first point (in degrees)
 ** \param[in]  lat1   Latitude of the first point (in degrees)
 ** \param[in]  long2  Longitude of the second point (in degrees)
 ** \param[in]  lat2   Latitude of the second point (in degrees)
 ** \param[in]  long3  Longitude of the third point (in degrees)
 ** \param[in]  lat3   Latitude of the third point (in degrees)
 **
 *****************************************************************************/
GEOSLIB_API double ut_geodetic_triangle_surface(double long1,
                                                double lat1,
                                                double long2,
                                                double lat2,
                                                double long3,
                                                double lat3)
{
  double a, b, c, A, B, C, surface;

  ut_geodetic_angles(long1, lat1, long2, lat2, long3, lat3, &a, &b, &c, &A, &B,
                     &C);
  surface = (A + B + C - GV_PI);
  return (surface);
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
GEOSLIB_API double ut_distance(int ndim, double *tab1, double *tab2)
{
  double distance, distang, R, v1, v2, delta;
  int flag_sphere;

  distance = 0.;
  variety_query(&flag_sphere);

  if (flag_sphere)
  {
    /* Case of the spherical coordinates */
    /* Longitude = 1st coord; Latitude = 2nd coord (in degrees) */

    variety_get_characteristics(&R);
    distang = ut_geodetic_angular_distance(tab1[0], tab1[1], tab2[0], tab2[1]);
    distance = R * distang;
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
GEOSLIB_API void ut_distance_allocated(int ndim, double **tab1, double **tab2)
{
  if (DISTANCE_NDIM < ndim)
  {
    DISTANCE_TAB1 = (double *) realloc((char *) DISTANCE_TAB1,
                                       sizeof(double) * ndim);
    DISTANCE_TAB2 = (double *) realloc((char *) DISTANCE_TAB2,
                                       sizeof(double) * ndim);
    DISTANCE_NDIM = ndim;
  }
  *tab1 = DISTANCE_TAB1;
  *tab2 = DISTANCE_TAB2;
  return;
}

/****************************************************************************/
/*!
 **  Calculate the intersection between two segments
 **
 ** \return 0 there is an intersection; 1 if there is no intersection
 **
 ** \param[in]  xd1,yd1     Starting point for the first segment
 ** \param[in]  xe1,ye1     Ending point for the first segment
 ** \param[in]  xd2,yd2     Starting point for the second segment
 ** \param[in]  xe2,ye2     Ending point for the second segment
 **
 ** \param[out]   xint,yint  Coordinates of the intersection
 **
 *****************************************************************************/
GEOSLIB_API int segment_intersect(double xd1,
                                  double yd1,
                                  double xe1,
                                  double ye1,
                                  double xd2,
                                  double yd2,
                                  double xe2,
                                  double ye2,
                                  double *xint,
                                  double *yint)
{
  double a1, a2, b1, b2, x, y, x1m, x1M, x2m, x2M, testval;

  /* Preliminary check */

  b1 = ye1 - yd1;
  b2 = ye2 - yd2;

  /* Case of two horizontal segments */

  if (ABS(b1) < EPSILON10 && ABS(b2) < EPSILON10)
  {
    if (ABS(ye1 - ye2) > EPSILON10) return (1);
    x1m = MIN(xd1, xe1);
    x1M = MAX(xd1, xe1);
    x2m = MIN(xd2, xe2);
    x2M = MAX(xd2, xe2);
    if (x1m > x2M || x2m > x1M) return (1);
    (*xint) = MAX(x1m, x2m);
    (*yint) = ye1;
    return (0);
  }

  /* Case of the horizontal first segment */

  if (ABS(b1) < EPSILON10)
  {
    y = ye1;
    x = xe2 + (y - ye2) * (xe2 - xd2) / b2;
    if ((x - xd1) * (x - xe1) > 0) return (1);
    if ((y - yd1) * (y - ye1) > 0) return (1);
    if ((x - xd2) * (x - xe2) > 0) return (1);
    if ((y - yd2) * (y - ye2) > 0) return (1);
    (*xint) = x;
    (*yint) = y;
    return (0);
  }

  /* Case of horizontal second segment */

  if (ABS(b2) < EPSILON10)
  {
    y = ye2;
    x = xe1 + (y - ye1) * (xe1 - xd1) / b1;
    if ((x - xd1) * (x - xe1) > 0) return (1);
    if ((y - yd1) * (y - ye1) > 0) return (1);
    if ((x - xd2) * (x - xe2) > 0) return (1);
    if ((y - yd2) * (y - ye2) > 0) return (1);
    (*xint) = x;
    (*yint) = y;
    return (0);
  }

  /* This operation is safe as end-point ordinates cannot be equal */

  a1 = (xe1 - xd1) / b1;
  a2 = (xe2 - xd2) / b2;

  /* Skip the case of parallel fractures */

  if (ABS(a1 - a2) < EPSILON10) return (1);
  y = (xd2 - xd1 + a1 * yd1 - a2 * yd2) / (a1 - a2);

  /* Discard intersection if located outside the segment */

  if (ABS(b1) > 0)
  {
    testval = (y - yd1) * (y - ye1);
    if (testval > 0) return (1);
  }
  if (ABS(b2) > 0)
  {
    testval = (y - yd2) * (y - ye2);
    if (testval > 0) return (1);
  }

  /* Update the endpoint in case of intersection */

  (*xint) = xd1 + a1 * (y - yd1);
  (*yint) = y;
  return (0);
}

/****************************************************************************/
/*!
 **  Evaluate the number of coefficients necessary to evaluate a function
 **  (at a sample location) at a given approximation
 **
 ** \return Minimum number of necessary coefficients
 **
 ** \param[in]  func       Function to be approximated
 ** \param[in]  cheb_elem  Cheb_Elem structure
 ** \param[in]  x          Sampling value
 ** \param[in]  nblin      Number of terms in the polynomial expansion
 ** \param[in]  blin       Array of coefficients for polynomial expansion
 **
 *****************************************************************************/
GEOSLIB_API int ut_chebychev_count(double (*func)(double,
                                                  double,
                                                  int,
                                                  double *),
                                   Cheb_Elem *cheb_elem,
                                   double x,
                                   int nblin,
                                   double *blin)
{
  double *coeffs, y, y0, T1, Tx, Tm1, Tm2, power, a, b, tol;
  int ncmax;

  // Initializations

  power = cheb_elem->power;
  a = cheb_elem->a;
  b = cheb_elem->b;
  tol = cheb_elem->tol;
  ncmax = cheb_elem->ncmax;
  coeffs = cheb_elem->coeffs;

  /* Get the true value */

  y0 = func(x, power, nblin, blin);

  /* Calculate the approximate value until tolerance is reached */

  T1 = 2 * (x - a) / (b - a) - 1.;
  y = coeffs[0] + coeffs[1] * T1;
  if (ABS(y * y - y0 * y0) / (y * y) < tol) return (2);
  Tm1 = T1;
  Tm2 = 1.;
  for (int i = 2; i < ncmax; i++)
  {
    Tx = 2. * T1 * Tm1 - Tm2;
    y += coeffs[i] * Tx;
    if (ABS(y * y - y0 * y0) / (y * y) < tol) return (i + 1);
    Tm2 = Tm1;
    Tm1 = Tx;
  }
  return (ncmax);
}

/*****************************************************************************/
/*!
 **  Calculates the coefficients of the Chebychev polynomial which is an
 **  approximation of a given function
 **
 **  \return  Error return code
 **
 ** \param[in]  func      Function to be approximated
 ** \param[in]  cheb_elem Cheb_Elem structure
 ** \param[in]  nblin     Number of terms in the polynomial expansion
 ** \param[in]  blin      Array of coefficients for polynomial expansion
 **
 *****************************************************************************/
GEOSLIB_API int ut_chebychev_coeffs(double (*func)(double,
                                                   double,
                                                   int,
                                                   double *),
                                    Cheb_Elem *cheb_elem,
                                    int nblin,
                                    double *blin)
{
  double *coeffs, *x1, *y1, *x2, *y2;
  double  minsubdiv, theta, ct, val1, val2, coeff, power, a, b;
  int     n, ncmax, error;

  /* Initializations */

  error = 1;
  power = cheb_elem->power;
  ncmax = cheb_elem->ncmax;
  a = cheb_elem->a;
  b = cheb_elem->b;
  coeffs = cheb_elem->coeffs;
  x1 = y1 = x2 = y2 = nullptr;

  minsubdiv = pow(2., 20.);
  if (minsubdiv >= (ncmax + 1) / 2)
    n = static_cast<int> (minsubdiv);
  else
    n = static_cast<int> (ceil((double) (ncmax + 1) / 2));

  /* Core allocation */

  x1 = (double *) mem_alloc(sizeof(double) * n, 0);
  if (x1 == nullptr) goto label_end;
  y1 = (double *) mem_alloc(sizeof(double) * n, 0);
  if (y1 == nullptr) goto label_end;
  x2 = (double *) mem_alloc(sizeof(double) * n, 0);
  if (x2 == nullptr) goto label_end;
  y2 = (double *) mem_alloc(sizeof(double) * n, 0);
  if (y2 == nullptr) goto label_end;

  /* Filling the arrays */

  for (int i = 0; i < n; i++)
  {
    theta = 2. * GV_PI * ((double) i) / ((double) n);
    ct = cos(theta / 2.);
    val1 = func(((b + a) + (b - a) * ct) / 2., power, nblin, blin);
    val2 = func(((b + a) - (b - a) * ct) / 2., power, nblin, blin);
    x1[i] = 0.5 * (val1 + val2);
    y1[i] = 0.;
    x2[i] = 0.5 * (val1 - val2) * cos(-theta / 2.);
    y2[i] = 0.5 * (val1 - val2) * sin(-theta / 2.);
  }

  /* Perform the FFT transform */

  if (fftn(1, &n, x1, y1,  1, 1.)) goto label_end;
  if (fftn(1, &n, x2, y2, -1, 1.)) goto label_end;

  /* Store the coefficients */

  coeff = 2. / (double) n;
  for (int i = 0; i < ncmax; i++)
    coeffs[i] = 0.;
  for (int i = 0; i < n; i++)
  {
    if (2 * i >= ncmax) break;
    coeffs[2 * i] = coeff * x1[i];
    if (2 * i + 1 >= ncmax) break;
    coeffs[2 * i + 1] = coeff * x2[i];
  }
  coeffs[0] /= 2.;

  /* Set the error return code */

  error = 0;

  label_end: x1 = (double *) mem_free((char * ) x1);
  y1 = (double *) mem_free((char * ) y1);
  x2 = (double *) mem_free((char * ) x2);
  y2 = (double *) mem_free((char * ) y2);
  return (error);
}

/****************************************************************************/
/*!
 **  Return the index of a sample when calculated from mirroring within
 **  an array whose indices vary between 0 and nx-1
 **
 ** \return Rank of the restrained cell
 **
 ** \param[in]  nx        Number of cells
 ** \param[in]  ix        Rank of the cell to be restrained
 **
 *****************************************************************************/
GEOSLIB_API int get_mirror_sample(int nx, int ix)
{
  int nmax;

  nmax = nx - 1;
  while (!(ix >= 0 && ix < nx))
  {
    if (ix < 0)
      ix = -ix;
    else if (ix > nmax) ix = 2 * nmax - ix;
  }
  return (ix);
}

/****************************************************************************/
/*!
 **  Rotation of a Direction in 3-D
 **
 ** \param[in]  ct,st     Cosine and Sine of the rotation angle
 ** \param[in]  a         Random direction
 ** \param[in,out] codir  Direction to be rotated
 **
 *****************************************************************************/
GEOSLIB_API void ut_rotation_direction(double ct,
                                       double st,
                                       double *a,
                                       double *codir)
{
  double rd, b[3], c[3], p[3];

  rd = 0.;
  for (int k = 0; k < 3; k++)
    rd += codir[k] * a[k];
  for (int k = 0; k < 3; k++)
    p[k] = rd * a[k];
  for (int k = 0; k < 3; k++)
    b[k] = codir[k] - p[k];

  rd = 0.;
  for (int k = 0; k < 3; k++)
    rd += b[k] * b[k];

  rd = sqrt(rd);
  for (int k = 0; k < 3; k++)
    b[k] /= rd;

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

  for (int k = 0; k < 3; k++)
    codir[k] = p[k] + rd * (ct * b[k] + st * c[k]);
}

/****************************************************************************/
/*!
 **  Generate a random rotation axis
 **
 ** \param[out]  ct       Cosine of the random rotation angle
 ** \param[out]  st       Sine of the random rotation angle
 ** \param[out]  a        Random direction vector
 **
 *****************************************************************************/
static void st_init_rotation(double *ct, double *st, double *a)
{
  double rd, theta;

  for (int k = 0; k < 3; k++)
    a[k] = law_gaussian();
  rd = 0.;
  for (int k = 0; k < 3; k++)
    rd += a[k] * a[k];
  rd = sqrt(rd);
  for (int k = 0; k < 3; k++)
    a[k] /= rd;

  theta = 2. * GV_PI * law_uniform(0., 1.);
  *ct = cos(theta);
  *st = sin(theta);
}

/****************************************************************************/
/*!
 **  Shuffle an array (by line)
 **
 ** \param[in]  nrow      Number of rows
 ** \param[in]  ncol      Number of columns
 ** \param[in,out] tab    Array to be suffled
 **
 *****************************************************************************/
GEOSLIB_API void ut_shuffle_array(int nrow, int ncol, double *tab)
{
  double *newtab, *rrank;
  int *irank, jrow;

  /* Core allocation */

  newtab = (double *) mem_alloc(sizeof(double) * nrow * ncol, 1);
  rrank = (double *) mem_alloc(sizeof(double) * nrow, 1);
  irank = (int *) mem_alloc(sizeof(int) * nrow, 1);

  /* Draw the permutation array */

  for (int i = 0; i < nrow; i++)
  {
    irank[i] = i;
    rrank[i] = law_uniform(0., 1.);
  }
  ut_sort_double(0, nrow, irank, rrank);

  /* Permutation from 'tab' into 'newtab' */

  for (int irow = 0; irow < nrow; irow++)
  {
    jrow = irank[irow];
    for (int icol = 0; icol < ncol; icol++)
      newtab[ncol * jrow + icol] = tab[ncol * irow + icol];
  }

  /* Restore in original array */

  for (int i = 0; i < nrow * ncol; i++)
    tab[i] = newtab[i];

  /* Core deallocation */

  irank = (int *) mem_free((char * ) irank);
  rrank = (double *) mem_free((char * ) rrank);
  newtab = (double *) mem_free((char * ) newtab);
}

/****************************************************************************/
/*!
 **  Generate a Van Der Corput list of points in R^3
 **
 ** \param[in]  n         Number of points
 ** \param[in]  flag_sym  Duplicate the samples by symmetry
 ** \param[in]  flag_rot  Perform a random rotation
 **
 ** \param[out] ntri_arg  Number of points
 ** \param[out] coor_arg  Array of point coordinates
 **                       (Dimension: 3*ntri)
 **
 *****************************************************************************/
GEOSLIB_API void ut_vandercorput(int n,
                                 int flag_sym,
                                 int flag_rot,
                                 int *ntri_arg,
                                 double **coor_arg)
{
  int i, j, ri, nb, ntri;
  double *coord, base, u, v, ct, st, a[3];

  /* Core allocation */

  ntri = 2 * n;
  coord = (double *) mem_alloc(sizeof(double) * 3 * ntri, 1);

  /* Processing */

  nb = 0;
  for (i = 0; i < n; i++)
  {

    // Binary decomposition
    j = i;
    u = 0;
    base = 2.;

    while (j)
    {
      ri = j % 2;
      u += ri / base;
      base *= 2;
      j = j / 2;
    }

    // Ternary decomposition
    j = i;
    v = 0;
    base = 3;

    while (j)
    {
      ri = j % 3;
      v += ri / base;
      base *= 3;
      j = j / 3;
    }

    COORD(0,nb) = cos(2. * GV_PI * u) * sqrt(1 - v * v);
    COORD(1,nb) = sin(2. * GV_PI * u) * sqrt(1 - v * v);
    COORD(2,nb) = v;
    nb++;

    if (flag_sym)
    {
      COORD(0,nb) = -cos(2. * GV_PI * u) * sqrt(1 - v * v);
      COORD(1,nb) = -sin(2. * GV_PI * u) * sqrt(1 - v * v);
      COORD(2,nb) = -v;
      nb++;
    }
  }

  /* Random rotation */

  if (flag_rot)
  {
    st_init_rotation(&ct, &st, a);
    for (i = 0; i < ntri; i++)
      ut_rotation_direction(ct, st, a, &coord[3 * i]);
  }

  /* Shuffle the samples in order to avoid three first colinear samples */

  ut_shuffle_array(ntri, 3, coord);

  /* Returning arguments */

  *ntri_arg = 2 * n;
  *coor_arg = coord;

  return;
}

/****************************************************************************/
/*!
 **  Set of auxiliary functions
 **
 *****************************************************************************/
/* normalize a vector of non-zero length */
static void st_normalize(double v[3])
{
  double d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] /= d;
  v[1] /= d;
  v[2] /= d;
}

/* Add the triangle */
static void st_addTriangle(double v1[3],
                           double v2[3],
                           double v3[3],
                           Reg_Coor *R_coor)
{
  int n;

  n = R_coor->ntri;

  R_coor->coor = (double *) mem_realloc((char * ) R_coor->coor,
                                        sizeof(double) * 3 * (n + 3), 1);

  for (int i = 0; i < 3; i++)
    RCOORD(i,n) = v1[i];
  for (int i = 0; i < 3; i++)
    RCOORD(i,n+1) = v2[i];
  for (int i = 0; i < 3; i++)
    RCOORD(i,n+2) = v3[i];
  R_coor->ntri += 3;
}

/* recursively subdivide face `depth' times */
void st_subdivide(double v1[3],
                  double v2[3],
                  double v3[3],
                  int depth,
                  Reg_Coor *R_coor)
{
  double v12[3], v23[3], v31[3];

  if (depth == 0)
  {
    st_addTriangle(v1, v2, v3, R_coor);
    return;
  }

  /* calculate midpoints of each side */
  for (int i = 0; i < 3; i++)
  {
    v12[i] = (v1[i] + v2[i]) / 2.0;
    v23[i] = (v2[i] + v3[i]) / 2.0;
    v31[i] = (v3[i] + v1[i]) / 2.0;
  }

  /* extrude midpoints to lie on unit sphere */
  st_normalize(v12);
  st_normalize(v23);
  st_normalize(v31);

  /* recursively subdivide new triangles */
  st_subdivide(v1, v12, v31, depth - 1, R_coor);
  st_subdivide(v2, v23, v12, depth - 1, R_coor);
  st_subdivide(v3, v31, v23, depth - 1, R_coor);
  st_subdivide(v12, v23, v31, depth - 1, R_coor);
}

static int st_already_present(Reg_Coor *R_coor, int i0, int ntri, double *coord)
{
  int found;
  static double eps = 1.e-03;

  if (ntri <= 0) return (0);
  for (int itri = 0; itri < ntri; itri++)
  {
    for (int k = found = 0; k < 3; k++)
      if (ABS(COORD(k,itri) - RCOORD(k,i0)) > eps) found = k + 1;
    if (found == 0) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Generate regular Icosahedron discretization
 **
 ** \return Error return code
 **
 ** \param[in]  n         Number of discretization steps
 ** \param[in]  flag_rot  Perform a random rotation
 **
 ** \param[out] ntri_arg  Number of points
 ** \param[out] coor_arg  Array of point coordinates
 **                       (Dimension: 3*ntri)
 **
 ** \remarks As random number are used in this function, a specific seed
 ** \remarks is fixed here
 **
 *****************************************************************************/
GEOSLIB_API int ut_icosphere(int n,
                             int flag_rot,
                             int *ntri_arg,
                             double **coor_arg)
{
  Reg_Coor R_coor;
  double *coord, ct, st, a[3];
  int ntri, seed_memo;
#define X 0.525731112119133696
#define Z 0.850650808352039932

  /* Get the current sed */

  seed_memo = law_get_random_seed();
  law_set_random_seed(43241);

  /* vertex data array */
  static double vdata[12][3] = { { -X, 0.0, Z },
                                 { X, 0.0, Z },
                                 { -X, 0.0, -Z },
                                 { X, 0.0, -Z },
                                 { 0.0, Z, X },
                                 { 0.0, Z, -X },
                                 { 0.0, -Z, X },
                                 { 0.0, -Z, -X },
                                 { Z, X, 0.0 },
                                 { -Z, X, 0.0 },
                                 { Z, -X, 0.0 },
                                 { -Z, -X, 0.0 } };

  /* triangle indices */
  static int tindices[20][3] = { { 1, 4, 0 }, { 4, 9, 0 }, { 4, 5, 9 }, { 8,
                                                                          5,
                                                                          4 },
                                 { 1, 8, 4 }, { 1, 10, 8 }, { 10, 3, 8 }, { 8,
                                                                            3,
                                                                            5 },
                                 { 3, 2, 5 }, { 3, 7, 2 }, { 3, 10, 7 }, { 10,
                                                                           6,
                                                                           7 },
                                 { 6, 11, 7 }, { 6, 0, 11 }, { 6, 1, 0 }, { 10,
                                                                            1,
                                                                            6 },
                                 { 11, 0, 9 }, { 2, 11, 9 }, { 5, 2, 9 }, { 11,
                                                                            2,
                                                                            7 } };

  if (n > 10)
  {
    messerr("The Regular Sphere discretization is limited to degree 10");
    law_set_random_seed(seed_memo);
    return (1);
  }
  R_coor.ntri = 0;
  R_coor.coor = nullptr;

  /* Subdivide the initial icosahedron */

  for (int i = 0; i < 20; i++)
  {
    st_subdivide(&vdata[tindices[i][0]][0], &vdata[tindices[i][1]][0],
                 &vdata[tindices[i][2]][0], n, &R_coor);
  }

  /* Suppress repeated triangle vertices */

  ntri = 0;
  coord = (double *) mem_alloc(sizeof(double) * 3 * R_coor.ntri, 1);
  for (int i = 0; i < R_coor.ntri; i++)
  {
    if (st_already_present(&R_coor, i, ntri, coord)) continue;
    for (int k = 0; k < 3; k++)
      COORD(k,ntri) = R_coor.coor[3 * i + k];
    ntri++;
  }

  /* Final resize */

  coord = (double *) mem_realloc((char * ) coord, sizeof(double) * 3 * ntri, 1);

  /* Random rotation */

  if (flag_rot)
  {
    st_init_rotation(&ct, &st, a);
    for (int i = 0; i < ntri; i++)
      ut_rotation_direction(ct, st, a, &coord[3 * i]);
  }

  /* Shuffle the samples in order to avoid three first colinear samples */

  ut_shuffle_array(ntri, 3, coord);

  /* Returning arguments */

  *ntri_arg = ntri;
  *coor_arg = coord;

  /* Free the Reg_Coor structure */

  R_coor.coor = (double *) mem_free((char * ) R_coor.coor);
  law_set_random_seed(seed_memo);
  return (0);
}

/****************************************************************************/
/*!
 **   Convert RGB into numeric
 **
 ** \param[in]  red     Red index
 ** \param[in]  green   Green index
 ** \param[in]  blue    Blue index
 **
 ** \param[out] c       Numeric value
 **
 *****************************************************************************/
GEOSLIB_API void rgb2num(int red,
                         int green,
                         int blue,
                         int /*a*/,
                         unsigned char *c)
{
  double value;

  value = (double) (red + green + blue) / 3.;

  if (value < 0.) value = 0.;
  if (value > 255.) value = 255.;
  *c = (unsigned char) value;

  return;
}

/****************************************************************************/
/*!
 **  Convert numeric to RGB
 **
 ** \param[in]   value   Input value
 **
 ** \param[out]  r     Red index
 ** \param[out]  g     Green index
 ** \param[out]  b     Blue index
 ** \param[out]  a     Transparency index
 **
 *****************************************************************************/
GEOSLIB_API void num2rgb(unsigned char value, int *r, int *g, int *b, int *a)
{
  *r = static_cast<int> ((value >> 24) & 0xff);
  *g = static_cast<int> ((value >> 16) & 0xff);
  *b = static_cast<int> ((value >>  8) & 0xff);
  *a = static_cast<int> ((value)       & 0xff);
}

/*****************************************************************************/
/*!
 **  Check if Legendre external functions have been defined properly
 **
 ** \return  1 if Legendre have been defined; 0 otherwise
 **
 ** \remarks This function returns a message is Legendre functions have not
 ** \remarks been defined
 **
 *****************************************************************************/
GEOSLIB_API int ut_is_legendre_defined(void)
{
  if (LEGENDRE_PL == nullptr)
  {
    messerr("You must define function 'legendre_Pl' beforehand");
    messerr("using the function 'define_legendre'");
    return (0);
  }
  if (LEGENDRE_SPHPLM == nullptr)
  {
    messerr("You must define function 'legendre_sphPlm' beforehand");
    messerr("using the function 'define_legendre'");
    return (0);
  }
  return (1);
}

/*****************************************************************************/
/*!
 **  Returns the Associated Legendre Function: legendre_Pl
 **
 ** \param[in]  flag_norm 1 for normalized and 0 otherwise
 ** \param[in]  n           Degree
 ** \param[in]  v           Value
 **
 *****************************************************************************/
GEOSLIB_API double ut_legendre(int flag_norm, int n, double v)
{
  int renard = -1;
  double res1 = 0.;
  double res2 = 0.;
  //double res3 = 0.;

  // TODO: Waiting for validation by Lantuejoul
  if (renard <= 0)
  {
    res1 = LEGENDRE_PL(n, v);
  }
//  if (renard == -1)
//  {
//    res3 = std::tr1::legendre(n, v);
//  }
  if (renard >= 0)
  {
    //  res2 = std::legendre(n,v);
    res2 = boost::math::legendre_p<double>(n, v);
  }

  if (renard == 0)
  {
    double diff = ABS(res1 + res2);
    if (diff > EPSILON5) diff = 100. * ABS(res1 - res2) / diff;
    if (diff > 5)
      messerr("---> Legendre n=%d v=%lf res1=%lf res2=%lf", n, v, res1, res2);
  }

  double result = res1;
  if (flag_norm)
  {
    double norme = sqrt((2. * ((double) n) + 1.) / 2.);
    result *= norme;
  }
  return (result);
}

/*****************************************************************************/
/*!
 **  Returns the Legendre Function legendre_Sphplm(n,k0,v) normalized
 **
 ** \param[in]  flag_norm 1 for normalized and 0 otherwise
 ** \param[in]  n           Degree
 ** \param[in]  k0          Order (ABS(k0) <= n)
 ** \param[in]  theta       Theta angle in radian
 **
 *****************************************************************************/
GEOSLIB_API double ut_flegendre(int flag_norm, int n, int k0, double theta)
{
  int k, flag_negative;
  int renard = -1;

  if (k0 < 0)
  {
    k = -k0;
    flag_negative = 1;
  }
  else
  {
    k = k0;
    flag_negative = 0;
  }

  // TODO: Waiting for the validation by Lantuejoul

  double res1 = 0.;
  double res2 = 0.;
  //double res3 = 0.;
  if (renard <= 0)
  {
    double v = cos(theta);
    res1 = LEGENDRE_SPHPLM(n, k, v);
    if (flag_negative && k % 2 == 1) res1 = -res1;
  }
//  if (renard == -1)
//  {
//    res3 = std::tr1::sph_legendre(n, k, theta);
//  }
//
  if (renard >= 0)
  {
    std::complex<double>
    resbis = boost::math::spherical_harmonic<double, double>(n, k, theta, 0.);
    res2 = resbis.real();
  }
  if (renard == 0)
  {
    double diff = ABS(res1 + res2);
    if (diff > EPSILON5) diff = 100. * ABS(res1 - res2) / diff;
    if (diff > 5)
      messerr(
          "---> Sph-Legendre n=%d k0=%d theta=%lf res1=%lf res2=%lf",
          n, k0, theta, res1, res2);
  }
  double result = res1;

  if (flag_norm)
  {
    double norme = 1. / sqrt(2 * GV_PI);
    result /= norme;
  }
  return (result);
}

/*****************************************************************************/
/*!
 **  Define the Legendre functions
 **
 ** \param[in]  legendre_sphPlm
 ** \param[in]  legendre_Pl
 **
 *****************************************************************************/
GEOSLIB_API void define_legendre(double (*legendre_sphPlm)(int, int, double),
                                 double (*legendre_Pl)(int, double))
{
  LEGENDRE_SPHPLM = legendre_sphPlm;
  LEGENDRE_PL = legendre_Pl;
}

/*****************************************************************************/
/*!
 **  Calculates the nbpoly log-factorial coefficients
 **
 ** \param[in]  nbpoly  Number of terms
 **
 ** \param[out] factor  logarithm of factorials
 **
 *****************************************************************************/
GEOSLIB_API void ut_log_factorial(int nbpoly, double *factor)
{
  int i;

  factor[0] = 0;
  for (i = 1; i < nbpoly; i++)
    factor[i] = factor[i - 1] + log((double) (i + 1));
}

/*****************************************************************************/
/*!
 **  Calculates the factorial coefficient
 **
 ** \return Returned value
 **
 ** \param[in]  k     Value
 **
 *****************************************************************************/
GEOSLIB_API double ut_factorial(int k)
{
  double val;

  val = 1;
  for (int i = 1; i <= k; i++)
    val *= (double) i;
  return (val);
}

/*****************************************************************************/
/*!
 **  Translates from degree to radian
 **
 ** \param[in]  angle  Angle in degrees
 **
 *****************************************************************************/
GEOSLIB_API double ut_deg2rad(double angle)
{
  return (angle * GV_PI / 180.);
}

/*****************************************************************************/
/*!
 **  Translates from radian to degree
 **
 ** \param[in]  angle  Angle in radian
 **
 *****************************************************************************/
GEOSLIB_API double ut_rad2deg(double angle)
{
  return (angle * 180. / GV_PI);
}

/****************************************************************************/
/*!
 **  Is a point inside a spherical triangle
 **
 ** \return 1 if the point belongs to the spherical triangle; 0 otherwise
 **
 ** \param[in]  coor    Coordinates of the target point (long,lat)
 ** \param[in]  surface Surface of the spherical triangle
 ** \param[in]  pts1    Coordinates of the first point of the triangle
 ** \param[in]  pts2    Coordinates of the second point of the triangle
 ** \param[in]  pts3    Coordinates of the third point of the triangle
 **
 ** \param[out] wgts    Array of weights
 **
 *****************************************************************************/
GEOSLIB_API int is_in_spherical_triangle(double *coor,
                                         double surface,
                                         double *pts1,
                                         double *pts2,
                                         double *pts3,
                                         double *wgts)
{
  double total, s[3], eps;

  eps = 1.e-6;
  total = 0.;
  s[0] = ut_geodetic_triangle_surface(coor[0], coor[1], pts2[0], pts2[1],
                                      pts3[0], pts3[1]);
  total += s[0];
  if (total > surface + eps) return (0);
  s[1] = ut_geodetic_triangle_surface(pts1[0], pts1[1], coor[0], coor[1],
                                      pts3[0], pts3[1]);
  total += s[1];
  if (total > surface + eps) return (0);
  s[2] = ut_geodetic_triangle_surface(pts1[0], pts1[1], pts2[0], pts2[1],
                                      coor[0], coor[1]);
  total += s[2];
  if (ABS(total - surface) > eps) return (0);
  for (int i = 0; i < 3; i++)
    wgts[i] = s[i] / total;
  return (1);
}

/****************************************************************************/
/*!
 ** Create a std::vector<double> for storing an array of double
 **
 ** \return The std::vector<double>
 **
 ** \param[in]  ntab      Number of samples
 ** \param[in]  rtab      Array of double values to be loaded
 **
 *****************************************************************************/
GEOSLIB_API std::vector<double> util_set_array_double(int ntab,
                                                      const double *rtab)
{
  if (debug_query("interface")) message("util_set_array_double\n");
  if (ntab <= 0 || rtab == nullptr) return std::vector<double>();
  std::vector<double> rettab(ntab);
  if (rettab.empty()) return rettab;

  for (int i = 0; i < ntab; i++)
    rettab[i] = rtab[i];

  return rettab;
}

/****************************************************************************/
/*!
 ** Create a std::vector<double> for storing an array of integer
 **
 ** \return  The std::vector<int>
 **
 ** \param[in]  ntab      Number of samples
 ** \param[in]  itab      Array of integer values to be loaded
 **
 *****************************************************************************/
GEOSLIB_API std::vector<int> util_set_array_integer(int ntab, const int *itab)
{
  if (debug_query("interface")) message("util_set_array_integer\n");
  std::vector<int> rettab(ntab);
  if (ntab <= 0 || itab == nullptr) return rettab;
  for (int i = 0; i < ntab; i++)
    rettab[i] = itab[i];
  return rettab;
}

/****************************************************************************/
/*!
 ** Create a std::vector<std::string> for storing an array of chars
 **
 ** \return  The std::vector<std::string>
 **
 ** \param[in]  ntab      Number of samples
 ** \param[in]  names     Array of character values to be loaded
 **
 *****************************************************************************/
GEOSLIB_API std::vector<std::string> util_set_array_char(int ntab, char **names)
{
  if (debug_query("interface")) message("util_set_array_char\n");
  std::vector<std::string> rettab(ntab);
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
GEOSLIB_API void set_last_message(int mode, const char *string)
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
        free((char *) LAST_MESSAGE[i]);
        LAST_MESSAGE[i] = nullptr;
      }
      free((char *) LAST_MESSAGE);
      NB_LAST_MESSAGE = 0;
      break;

    case 1:                       // Add string to array of messages
      size = static_cast<int> (strlen(string));
      if (size <= 0) return;

      if (NB_LAST_MESSAGE <= 0)
        LAST_MESSAGE = (char **) malloc(sizeof(char *) * 1);
      else
        LAST_MESSAGE = (char **) realloc(
            (char *) LAST_MESSAGE, sizeof(char *) * (NB_LAST_MESSAGE + 1));
      LAST_MESSAGE[NB_LAST_MESSAGE] = address = (char *) malloc(size + 1);
      (void) gslStrcpy(address, string);
      address[size] = '\0';
      NB_LAST_MESSAGE++;
      break;

    case -1:                    // Concatenate
      size = static_cast<int> (strlen(string));
      if (size <= 0) return;

      if (NB_LAST_MESSAGE <= 0)
      {
        set_last_message(1, string);
        return;
      }

      sizaux = static_cast<int> (strlen(LAST_MESSAGE[NB_LAST_MESSAGE - 1]));
      LAST_MESSAGE[NB_LAST_MESSAGE - 1] = address = (char *) realloc(
          (char *) LAST_MESSAGE[NB_LAST_MESSAGE - 1], size + sizaux + 2);
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
GEOSLIB_API void print_last_message(void)
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
 **  Return all the combinations of k within n (local recursive routine)
 **
 ** \param[in]  v     Array of indices to be sorted
 ** \param[in]  start Rank of the starting index
 ** \param[in]  n     Total number of objects (>= 1)
 ** \param[in]  k     Starting sorting index
 ** \param[in]  maxk  Selected number of objects (>= 1)
 **
 ** \param[in,out] ncomb  Current number of combinations
 ** \param[in,out] comb   Current array of combinations
 **
 *****************************************************************************/
static void st_combinations(int *v,
                            int start,
                            int n,
                            int k,
                            int maxk,
                            int *ncomb,
                            int **comb)
{
  int i, nloc, ndeb, *cloc;

  nloc = *ncomb;
  cloc = *comb;

  /* k here counts through positions in the maxk-element v.
   * if k > maxk, then the v is complete and we can use it.
   */
  if (k > maxk)
  {
    /* insert code here to use combinations as you please */
    cloc = (int *) mem_realloc((char * ) cloc, sizeof(int) * maxk * (nloc + 1),
                               1);
    ndeb = nloc * maxk;
    for (i = 0; i < maxk; i++)
      cloc[ndeb + i] = v[i + 1];
    *ncomb = nloc + 1;
    *comb = cloc;
    return;
  }

  /* for this k'th element of the v, try all start..n
   * elements in that position
   */
  for (i = start; i <= n; i++)
  {
    v[k] = i;

    /* recursively generate combinations of integers from i+1..n
     */
    st_combinations(v, i + 1, n, k + 1, maxk, ncomb, comb);
  }
}

/****************************************************************************/
/*!
 **  Return all the combinations of k within n
 **
 ** \return Return all the combinations of 'k' objects amongst 'n'
 **
 ** \param[in]  n     Total number of objects (>1)
 ** \param[in]  maxk  Selected number of objects (1<=maxk<n)
 **
 ** \param[out] ncomb Number of combinations
 **
 ** \remarks The calling function must free the returned array.
 **
 *****************************************************************************/
GEOSLIB_API int *ut_combinations(int n, int maxk, int *ncomb)
{
  int *v, *comb;

  v = (int *) mem_alloc(sizeof(int) * n, 1);
  for (int i = 0; i < n; i++)
    v[i] = i;

  (*ncomb) = 0;
  comb = nullptr;
  st_combinations(v, 1, n, 1, maxk, ncomb, &comb);
  v = (int *) mem_free((char * ) v);
  return (comb);
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
GEOSLIB_API int *ut_split_into_two(int ncolor,
                                   int flag_half,
                                   int verbose,
                                   int *nposs)
{
  int p, nmax, ncomb, np, lec;
  int *mattab, *comb;

  /* Initializations */

  p = (flag_half) ? static_cast<int> (floor((double) ncolor / 2.)) : ncolor - 1;
  nmax = static_cast<int> (pow(2, ncolor));
  mattab = comb = nullptr;
  np = 0;

  /* Core allocation */

  mattab = (int *) mem_alloc(sizeof(int) * ncolor * nmax, 1);
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
  comb = (int *) mem_free((char * ) comb);

  /* Resize */

  mattab = (int *) mem_realloc((char * ) mattab, sizeof(int) * ncolor * np, 1);
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
 **  Is a point inside a spherical triangle
 **
 ** \return 1 if the point belongs to the spherical triangle; 0 otherwise
 **
 ** \param[in]  coor    Coordinates of the target point (long,lat)
 ** \param[in]  ptsa    Coordinates of the first point of the triangle
 ** \param[in]  ptsb    Coordinates of the second point of the triangle
 ** \param[in]  ptsc    Coordinates of the third point of the triangle
 **
 ** \param[out] wgts    Array of weights
 **
 *****************************************************************************/
GEOSLIB_API int is_in_spherical_triangle_optimized(double *coor,
                                                   double *ptsa,
                                                   double *ptsb,
                                                   double *ptsc,
                                                   double *wgts)
{
  double total, s[3], stot, eps;
  double A, B, C, AB, AC, BA, BC, CA, CB, OA, OB, OC;
  double dab, dbc, dac, d0a, d0b, d0c;
  double sinab, cosab, sinbc, cosbc, sinac, cosac;
  double sin0a, cos0a, cos0b, sin0b, cos0c, sin0c;

  eps = 1.e-6;
  total = 0.;

  dab = ut_geodetic_angular_distance(ptsa[0], ptsa[1], ptsb[0], ptsb[1]);
  dbc = ut_geodetic_angular_distance(ptsb[0], ptsb[1], ptsc[0], ptsc[1]);
  dac = ut_geodetic_angular_distance(ptsa[0], ptsa[1], ptsc[0], ptsc[1]);
  d0a = ut_geodetic_angular_distance(coor[0], coor[1], ptsa[0], ptsa[1]);
  d0b = ut_geodetic_angular_distance(coor[0], coor[1], ptsb[0], ptsb[1]);
  d0c = ut_geodetic_angular_distance(coor[0], coor[1], ptsc[0], ptsc[1]);

  sinab = sin(dab);
  cosab = cos(dab);
  sinbc = sin(dbc);
  cosbc = cos(dbc);
  sinac = sin(dac);
  cosac = cos(dac);
  sin0a = sin(d0a);
  cos0a = cos(d0a);
  sin0b = sin(d0b);
  cos0b = cos(d0b);
  sin0c = sin(d0c);
  cos0c = cos(d0c);

  A = st_convert_geodetic_angle(sinbc, cosbc, sinac, cosac, sinab, cosab);
  B = st_convert_geodetic_angle(sinac, cosac, sinab, cosab, sinbc, cosbc);
  C = st_convert_geodetic_angle(sinab, cosab, sinbc, cosbc, sinac, cosac);
  stot = (A + B + C - GV_PI);

  OA = st_convert_geodetic_angle(sinbc, cosbc, sin0c, cos0c, sin0b, cos0b);
  BA = st_convert_geodetic_angle(sin0c, cos0c, sin0b, cos0b, sinbc, cosbc);
  CA = st_convert_geodetic_angle(sin0b, cos0b, sinbc, cosbc, sin0c, cos0c);
  s[0] = (OA + BA + CA - GV_PI);
  total += s[0];
  if (total > stot + eps) return (0);

  AB = st_convert_geodetic_angle(sin0c, cos0c, sinac, cosac, sin0a, cos0a);
  OB = st_convert_geodetic_angle(sinac, cosac, sin0a, cos0a, sin0c, cos0c);
  CB = st_convert_geodetic_angle(sin0a, cos0a, sin0c, cos0c, sinac, cosac);
  s[1] = (AB + OB + CB - GV_PI);
  total += s[1];
  if (total > stot + eps) return (0);

  AC = st_convert_geodetic_angle(sin0b, cos0b, sin0a, cos0a, sinab, cosab);
  BC = st_convert_geodetic_angle(sin0a, cos0a, sinab, cosab, sin0b, cos0b);
  OC = st_convert_geodetic_angle(sinab, cosab, sin0b, cos0b, sin0a, cos0a);
  s[2] = (AC + BC + OC - GV_PI);
  total += s[2];
  if (ABS(total - stot) > eps) return (0);

  for (int i = 0; i < 3; i++)
    wgts[i] = s[i] / total;
  return (1);
}

/****************************************************************************/
/*!
 **  Decode the grid sorting order
 **
 ** \return Array describing the order
 **
 ** \param[in]  name    Name of the sorting string
 ** \param[in]  ndim    Space dimension
 ** \param[in]  nx      Array giving the number of cells per direction
 ** \param[in]  verbose Verbose flag
 **
 ** \remarks The value of order[i] gives the dimension of the space along
 ** \remarks which sorting takes place at rank "i". This value is positive
 ** \remarks for increasing order and negative for decreasing order.
 ** \remarks 'order' values start from 1 (for first space dimension)
 **
 ** \remarks The calling function must free the returned array
 **
 *****************************************************************************/
GEOSLIB_API int *ut_name_decode(const char *name,
                                int ndim,
                                int *nx,
                                int verbose)
{
  int *order, *ranks, num, orient, idim, error, a_order;
  char *p;

  // Initializations 

  error = 1;
  order = ranks = nullptr;

  // Core allocation

  order = (int *) mem_alloc(sizeof(int) * ndim, 0);
  if (order == nullptr) goto label_end;
  ranks = (int *) mem_alloc(sizeof(int) * ndim, 0);
  if (ranks == nullptr) goto label_end;
  for (int i = 0; i < ndim; i++)
    ranks[i] = 0;

  // Loop on the character string

  idim = 0;
  p = (char *) name;
  while (*p)
  {
    if ((*p) == '-' && (*(p + 1)) == 'x' && isdigit(*(p + 2)))
      orient = -1;
    else if ((*p) == '+' && (*(p + 1)) == 'x' && isdigit(*(p + 2)))
      orient = 1;
    else
      orient = 0;

    if (orient != 0)
    {

      // The string '+x*' or '-x*' has been encountered

      p += 2;
      num = (int) strtol(p, &p, 10);
      order[idim] = orient * num;
      if (num > ndim)
      {
        messerr("'order' refers to 'x%d' while space dimension is %d", num,
                ndim);
        goto label_end;
      }
      ranks[num - 1] = orient;
      idim++;
    }
    else
    {
      p++;
    }
  }

  // Check that all indices (below space dimension) have been specified

  for (int i = 0; i < ndim; i++)
  {
    if (ranks[i] != 0) continue;
    messerr("'x%d' is not mentionned in 'order'", i + 1);
    goto label_end;
  }

  // Optional printout

  if (verbose)
  {
    message("Decoding the sorting rule (%s) with nx = (", name);
    for (int i = 0; i < ndim; i++)
      message(" %d", nx[i]);
    message(" )\n");
    for (int i = 0; i < ndim; i++)
    {
      a_order = ABS(order[i]);
      message("%d - Dimension=%d - N%d=%d", i + 1, a_order, a_order,
              nx[a_order - 1]);
      if (order[i] > 0)
        message(" - Increasing\n");
      else
        message(" - Decreasing\n");
    }
  }

  // Set the error return code

  error = 0;

  label_end: ranks = (int *) mem_free((char * ) ranks);
  if (error) order = (int *) mem_free((char * ) order);
  return (order);
}

/****************************************************************************/
/*!
 ** Perform the recursion through the dimension
 **
 ** \param[in]  idim    Space rank
 ** \param[in]  verbose Verbose flag
 ** \param[in]  int_str Pointer to the internal structure
 **
 *****************************************************************************/
static void st_dimension_recursion(int idim, int verbose, void *int_str)
{
  Dim_Loop *dlp;
  int order, sdim, nval, ival, ndim;

  // Assignments

  dlp = (Dim_Loop *) int_str;
  ndim = dlp->ndim;

  if (idim < 0)
  {

    // We have reached the bottom of the pile, evaluate the absolute address

    ival = dlp->indg[ndim - 1];
    for (int jdim = ndim - 2; jdim >= 0; jdim--)
      ival = ival * dlp->nx[jdim] + dlp->indg[jdim];
    dlp->tab[dlp->curech++] = ival + 1;

    // Optional printout

    if (verbose)
    {
      message("node (");
      for (int jdim = 0; jdim < ndim; jdim++)
        message(" %d", dlp->indg[jdim] + 1);
      message(" ) -> %d\n", ival + 1);
    }
    return;
  }

  order = dlp->order[idim];
  sdim = ABS(order) - 1;
  nval = dlp->nx[sdim];

  // Loop

  for (int jy = 0; jy < nval; jy++)
  {
    dlp->indg[sdim] = (order < 0) ? nval - jy - 1 :
                                    jy;
    st_dimension_recursion(idim - 1, verbose, (void *) dlp);
  }
  return;
}

/****************************************************************************/
/*!
 **  Allocates and returns an array giving the ranks of the
 **  cells (sequentially according to user's order) coded
 **  with standard ranks (according to Geoslib internal order)
 **
 ** \return Array of indices
 **
 ** \param[in]  ndim    Space dimension
 ** \param[in]  nx      Array giving the number of cells per direction
 ** \param[in]  order   Array of grid orders
 ** \param[in]  verbose Verbose flag
 **
 *****************************************************************************/
GEOSLIB_API double *ut_rank_cells(int ndim, int *nx, int *order, int verbose)
{
  double *tab, *tab2;
  int *indg, *ind, error, ncell;
  Dim_Loop dlp;

  // Initializations

  error = 1;
  tab = tab2 = nullptr;
  indg = ind = nullptr;
  ncell = 1;
  for (int idim = 0; idim < ndim; idim++)
    ncell *= nx[idim];

  // Core allocation

  tab = (double *) mem_alloc(sizeof(double) * ncell, 0);
  if (tab == nullptr) goto label_end;
  indg = (int *) mem_alloc(sizeof(int) * ndim, 0);
  if (indg == nullptr) goto label_end;
  ind = (int *) mem_alloc(sizeof(int) * ncell, 0);
  if (ind == nullptr) goto label_end;
  tab2 = (double *) mem_alloc(sizeof(double) * ncell, 0);
  if (tab2 == nullptr) goto label_end;
  for (int i = 0; i < ndim; i++)
    indg[i] = 0;

  // Initialization of the recursion structure

  dlp.curech = 0;
  dlp.ndim = ndim;
  dlp.nx = nx;
  dlp.order = order;
  dlp.indg = indg;
  dlp.tab = tab;

  // Recursion

  st_dimension_recursion(ndim - 1, verbose, (void *) &dlp);

  // Invert order

  for (int i = 0; i < ncell; i++)
    ind[i] = i;
  ut_sort_double(0, ncell, ind, tab);
  for (int i = 0; i < ncell; i++)
    tab2[i] = tab[ind[i]];

  // Set the error returned code

  error = 0;

  label_end: if (error) tab = (double *) mem_free((char * ) tab);
  tab = (double *) mem_free((char * ) tab);
  ind = (int *) mem_free((char * ) ind);
  indg = (int *) mem_free((char * ) indg);
  return (tab2);
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
char *convert(const std::string & s)
{
  char *pc = new char[s.size() + 1];
  std::strcpy(pc, s.c_str());
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
GEOSLIB_API std::vector<char *> util_vs_to_vs(VectorString vs)
{
  std::vector<char *> vc;
  std::transform(vs.begin(), vs.end(), std::back_inserter(vc), convert);
  return vc;
}

/****************************************************************************/
/*!
 **   Convert angles to rotation matrices
 **
 ** \param[in]  ndim     Number of space dimensions
 ** \param[in]  ndir     Number of directions
 ** \param[in]  angles   Vector giving the angles characteristics (in degrees)
 **
 ** \param[out] codir    Vector of the direction (Dim: ndir * ndim)
 **
 ** \remarks If angles is not provided:
 ** \remarks - if ndir == ndim: return basic direction of space
 ** \remarks - if ndim < ndim: return 'ndir' directions regular in the 2-D
 **
 *****************************************************************************/
GEOSLIB_API void ut_angles_to_codir(int ndim,
                                    int ndir,
                                    const VectorDouble& angles,
                                    VectorDouble& codir)
{
  if (ndim <= 1) return;

  codir.resize(ndim * ndir);
  for (int i = 0; i < ndim * ndir; i++) codir[i] = 0.;

  if (angles.size() <= 0)
  {
    if (ndir == ndim)
    {
      int ecr = 0;
      for (int idir = 0; idir < ndir; idir++)
        for (int idim = 0; idim < ndim; idim++)
          codir[ecr++] = (idir == idim) ? 1. : 0.;
    }
    else
    {
      {
        for (int idir = 0; idir < ndir; idir++)
        {
          double angref = 180. * idir / ndir;
          codir[idir * ndim + 0] = cos(angref * GV_PI / 180.);
          codir[idir * ndim + 1] = sin(angref * GV_PI / 180.);
        }
      }
    }
  }
  else
  {
    if ((int) angles.size() == ndir)
    {
      for (int idir = 0; idir < ndir; idir++)
      {
        codir[idir * ndim + 0] = cos(angles[idir] * GV_PI / 180.);
        codir[idir * ndim + 1] = sin(angles[idir] * GV_PI / 180.);
      }
    }
  }
  return;
}

/****************************************************************************/
/*!
 **   Perform a search of a Pattern within a String
 **   The input String can contain regular expression
 **
 ** \return  1 if there is a match; 0 otherwise
 **
 ** \param[in]  string    String to be serached
 ** \param[in]  pattern   Search Pattern
 ** \param[in]  verbose   Verbose flag
 **
 *****************************************************************************/
static int st_string_search(const String& string,
                            const String& pattern,
                            int verbose)
{
  int ok = 0;

  try
  {
    std::regex string_regex(pattern);
    ok = (std::regex_match(string, string_regex));
    if (verbose)
    {
      message("Searching '%s' in '%s' : ",pattern.c_str(),string.c_str());
      if (ok)
        message("OK\n");
      else
        message("Not found\n");
    }
  }
  catch (std::regex_error& e)
  {
    std::cerr << "Invalid Regular Expression." << std::endl;
  }
  return ok;
}

/****************************************************************************/
/*!
 **   Perform a search of a Pattern within a VectorString
 **   and returns the list of matches.
 **   The input String can contain regular expression
 **
 ** \return  VectorInt of matching indices (starting from 1)
 **
 ** \param[in]  list_string     VectorString structure
 ** \param[in]  pattern         Search Pattern
 ** \param[in]  verbose         Verbose flag
 **
 *****************************************************************************/
GEOSLIB_API VectorInt util_string_search(const VectorString& list_string,
                                         const String& pattern,
                                         int verbose)
{
  VectorInt ranks;
  int ns = static_cast<int> (list_string.size());
  for (int is = 0; is < ns; is++)
  {
    if (st_string_search(list_string[is], pattern, verbose))
      ranks.push_back(is+1);
  }
  return ranks;
}
