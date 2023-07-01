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
#include "geoslib_d.h"
#include "geoslib_old_f.h"

#include "Enum/EDbg.hpp"

#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/Law.hpp"

#include <cmath>
#include <map>

#define LSTACK    1000
#define MINI        10

static EDbg _debugOptions = EDbg::DB;

bool isInteger(double value, double eps)
{
  int iclose = getClosestInteger(value);
  if (ABS((double) iclose - value) > eps) return false;
  return true;
}

int getClosestInteger(double value)
{
  int iclose = (int) round(value);
  return iclose;
}

bool isMultiple(int nbig, int nsmall)
{
  double ratio;

  ratio = (double) nbig / (double) nsmall;
  return (isInteger(ratio));
}

bool isOdd(int number)
{
  int middle;

  middle = number / 2;
  if (number != 2 * middle)
    return true;
  else
    return false;
}

bool isEven(int number)
{
  int middle;

  middle = number / 2;
  if (number != 2 * middle)
    return false;
  else
    return true;
}

double getMin(double val1, double val2)
{
  if (FFFF(val1)) return (val2);
  if (FFFF(val2)) return (val1);
  return (MIN(val1, val2));
}

double getMax(double val1, double val2)
{
  if (FFFF(val1)) return (val2);
  if (FFFF(val2)) return (val1);
  return (MAX(val1, val2));
}

#ifndef SWIG

double getTEST()
{
  return TEST;
}

int getITEST()
{
  return ITEST;
}

/****************************************************************************/
/*!
 **  Checks if a double value is TEST
 **
 ** \return  1 if a TEST value is encountered; 0 otherwise
 **
 ** \param[in]  value Value to be tested
 **
 *****************************************************************************/
int FFFF(double value)
{
  int rep;

  rep = 0;
  if (std::isnan(value)) rep = 1; // TODO : what about std::isinf ?
  if (value > TEST_COMP) rep = 1;

  return (rep);
}

/****************************************************************************/
/*!
 **  Checks if an integer value is TEST
 **
 ** \return  1 if a ITEST value is encountered; 0 otherwise
 **
 ** \param[in]  value Value to be tested
 **
 *****************************************************************************/
int IFFFF(int value)
{
  int rep;

  rep = 0;
  if (value == ITEST) rep = 1;

  return (rep);
}

#endif //SWIG

/*****************************************************************************/
/*!
 **  Translates from degree to radian
 **
 ** \param[in]  angle  Angle in degrees
 **
 *****************************************************************************/
double ut_deg2rad(double angle)
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
double ut_rad2deg(double angle)
{
  return (angle * 180. / GV_PI);
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
void ut_sort_double(int safe, int nech, int *ind, double *value)
{
  static int LISTE_L[LSTACK];
  static int LISTE_R[LSTACK];
  int i, j, p, l, r, pstack, inddev, inddeu;
  double *tab, tablev, tableu;

  /* Initialization */

  inddev = inddeu = 0;
  if (safe)
  {
    tab = (double*) mem_alloc(sizeof(double) * nech, 1);
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

  if (safe) tab = (double*) mem_free((char* ) tab);
  return;
}

/****************************************************************************/
/*!
 **  Returns the statistics of an array in a StatResults structure
 **
 ** \param[in]  nech    Number of samples
 ** \param[in]  tab     Array of values
 ** \param[in]  sel     Array containing the Selection or NULL
 ** \param[in]  wgt     Array containing the Weights or NULL
 **
 ****************************************************************************/
StatResults ut_statistics(int nech, double *tab, double *sel, double *wgt)
{
  StatResults stats;

  /* Initializations */

  double tmin =  1.e30;
  double tmax = -1.e30;
  double num = 0.;
  double mm = 0.;
  double vv = 0.;
  int nval = 0;

  for (int i = 0; i < nech; i++)
  {
    if (sel != nullptr && sel[i] == 0.) continue;
    if (FFFF(tab[i])) continue;
    double weight = (wgt != nullptr && wgt[i] >= 0) ? wgt[i] : 1.;
    if (tab[i] < tmin) tmin = tab[i];
    if (tab[i] > tmax) tmax = tab[i];
    nval++;
    num += weight;
    mm  += weight * tab[i];
    vv  += weight * tab[i] * tab[i];
  }

  /* Returning arguments */

  stats.number = nech;
  stats.nvalid = nval;
  if (tmax < tmin || nval <= 0)
  {
    stats.mini = TEST;
    stats.maxi = TEST;
    stats.delta = TEST;
    stats.mean = TEST;
    stats.stdv = TEST;
  }
  else
  {
    stats.mini = tmin;
    stats.maxi = tmax;
    stats.delta = tmax - tmin;
    mm /= num;
    vv = vv / num - mm * mm;
    if (vv < 0.) vv = 0.;
    stats.mean = mm;
    stats.stdv = sqrt(vv);
  }
  return stats;
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
void ut_stats_mima_print(const char *title, int nech, double *tab, double *sel)
{
  StatResults stats = ut_statistics(nech, tab, sel);

  // Print the statistics out

  if (stats.nvalid <= 0)
    message("%s: NVal=%6d/%6d - Min=NA - Max=NA\n", title, stats.nvalid, nech);
  else
    message("%s: NVal=%6d/%6d - Min=%lf - Max=%lf\n", title, stats.nvalid, nech,
            stats.mini, stats.maxi);
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
void ut_facies_statistics(int nech,
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
    *nval = 0;
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
void ut_classify(int nech,
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
double ut_median(double *tab, int ntab)
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
 **  Compute combinations(n,k)
 **
 ** \return Return the number of combinations of 'k' objects amongst 'n'
 **
 ** \param[in]  n     Total number of objects (>= 1)
 ** \param[in]  k     Selected number of objects (>= 1)
 **
 *****************************************************************************/
double ut_cnp(int n, int k)
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
double* ut_pascal(int ndim)
{
  double *m;
#define M(j,i)            (m[(i) * ndim + (j)])

  /* Core allocation */

  m = (double*) mem_alloc(sizeof(double) * ndim * ndim, 0);
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
    cloc = (int*) mem_realloc((char* ) cloc, sizeof(int) * maxk * (nloc + 1),
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
int* ut_combinations(int n, int maxk, int *ncomb)
{
  int *v, *comb;

  v = (int*) mem_alloc(sizeof(int) * n, 1);
  for (int i = 0; i < n; i++)
    v[i] = i;

  (*ncomb) = 0;
  comb = nullptr;
  st_combinations(v, 1, n, 1, maxk, ncomb, &comb);
  v = (int*) mem_free((char* ) v);
  return (comb);
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
void ut_shuffle_array(int nrow, int ncol, double *tab)
{
  int jrow;

  /* Core allocation */

  VectorDouble newtab(nrow * ncol);
  VectorDouble rrank(nrow);
  VectorInt irank(nrow);

  /* Draw the permutation array */

  for (int i = 0; i < nrow; i++)
  {
    irank[i] = i;
    rrank[i] = law_uniform(0., 1.);
  }
  VH::arrangeInPlace(0, irank, rrank, true, nrow);

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
}

/**
 * Returns the list of absolute indices for the only active samples
 * A sample is active if its 'sel' value is equal to 1
 * @param sel Vector giving the status of all samples (Dimension: absolute)
 * @return
 */
VectorInt getListActiveToAbsolute(const VectorDouble& sel)
{
  int nech = (int) sel.size();
  VectorInt ranks;
  for (int iabs = 0; iabs < nech; iabs++)
  {
    if (sel[iabs]) ranks.push_back(iabs);
  }
  return ranks;
}

/**
 * Returns the map such that MAP[iabs] = iact.
 * A sample is active if its 'sel' value is equal to 1
 * @param sel Vector giving the status of all samples (Dimension: absolute)
 * @param verbose Verbose flag
 * @return The map (dimension: nrel)
 */
std::map<int, int> getMapAbsoluteToRelative(const VectorDouble& sel, bool verbose)
{
  std::map<int, int> map;
  int nabs = (int) sel.size();
  int ifirst = ITEST;
  int ilast  = ITEST;
  int irel   = 0;
  for (int iabs = 0; iabs < nabs; iabs++)
  {
    if (sel[iabs] == 0) continue;
    map[iabs] = irel++;

    if (IFFFF(ifirst)) ifirst = iabs;
    ilast = iabs;
  }

  // Optional control printout
  if (verbose)
  {
    message("Map Absolute to Relative\n");
    message("- Number of absolute positions = %d\n", nabs);
    message("- Number of active positions   = %d\n", irel);
    message("- Absolute address of the first active sample = %d\n",ifirst);
    message("- Absolute address of the last active sample  = %d\n",ilast);
  }
  return map;
}

/**
 * Returns the rank of the relative grid node from its absolute index using the Map
 * @param map  The <int,int> map
 * @param iabs Absolute rank of the grid node
 * @return Rank of the corresponding active (relative) grid node (or -1 is not found)
 */
int getRankMapAbsoluteToRelative(const std::map<int, int>& map, int iabs)
{
  if (map.empty()) return iabs;
  if (map.find(iabs) == map.end())
    return -1;
  else
    return map.find(iabs)->second;
}

int getRankMapRelativeToAbsolute(const std::map<int, int>& map, int irel)
{
  if (map.empty()) return irel;
  auto it = map.begin();
  std::advance(it, irel);
  return it->first;
}
