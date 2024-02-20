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
#include "Db/DbHelper.hpp"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Skin/Skin.hpp"
#include "Basic/Law.hpp"

#include "geoslib_old_f.h"

static DbGrid *DB_GRID_FILL;

#define R(i,j)              (R[(i) * n + (j)])

class LocalSkin: public ISkinFunctions
{
  /****************************************************************************/
  /*!
   **  Check if the cell is already filled
   **
   ** \return  1 if the cell (filled with facies) is already filled
   **
   ** \param[in]  ipos  Absolute grid index of the input grid node
   **
   *****************************************************************************/
  int isAlreadyFilled(int ipos) const override
  {
    if (!DB_GRID_FILL->getSelection(ipos)) return (0);
    int value = FFFF(DB_GRID_FILL->getLocVariable(ELoc::Z,ipos, 0)) ? 0 : 1;
    return (value);
  }
  /****************************************************************************/
  /*!
   **  Check if the cell can be filled with fluid
   **
   ** \return  1 if the cell (filled with facies) can be filled with Fluid
   **
   ** \param[in]  ipos  Absolute grid index of the input grid node
   **
   *****************************************************************************/
  int isToBeFilled(int ipos) const override

  {
    if (!DB_GRID_FILL->getSelection(ipos)) return (0);
    int value = FFFF(DB_GRID_FILL->getLocVariable(ELoc::Z,ipos, 0)) ? 1 : 0;
    return (value);
  }
};

/****************************************************************************/
/*!
 **  Check if a pair must be kept according to code criterion
 **
 ** \return  1 if the codes are not comparable
 **
 ** \param[in]  db1        First Db structure
 ** \param[in]  db2        Second Db structure
 ** \param[in]  iech       Rank of the first sample
 ** \param[in]  jech       Rank of the second sample
 ** \param[in]  opt_code   code selection option
 ** \li                     0 : no use of the code selection
 ** \li                     1 : codes must be close enough
 ** \li                     2 : codes must be different
 ** \param[in]  tolcode    Code tolerance
 **
 ** \remarks When used in variogram calculation, pairs are discarded then the
 ** \remarks resulting value is 1.
 **
 *****************************************************************************/
static int st_code_comparable(const Db *db1,
                              const Db *db2,
                              int iech,
                              int jech,
                              int opt_code,
                              int tolcode)
{
  double code1, code2;

  /* Dispatch */

  switch (opt_code)
  {
    case 0:
      break;

    case 1: /* Code must be close */
      code1 = db1->getLocVariable(ELoc::C,iech,0);
      code2 = db2->getLocVariable(ELoc::C,jech,0);
      if (ABS(code1 - code2) > tolcode) return (1);
      break;

    case 2: /* Code must be different */
      code1 = db1->getLocVariable(ELoc::C,iech,0);
      code2 = db2->getLocVariable(ELoc::C,jech,0);
      if (code1 == code2) return (1);
      break;
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Find the neighborhood of the current cell
 **
 ** \param[in]  ipos    Absolute grid index of the input grid node
 ** \param[in]  ndim    Space dimension
 ** \param[in]  radius  Radius of the neighborhood
 **
 ** \param[out] nech_loc Number of samples in the neighborhood
 ** \param[out] tabind   Index of the neighboring sample
 ** \param[out] tabval   Value of the neighboring sample
 **
 *****************************************************************************/
static void st_grid_fill_neigh(int ipos,
                               int ndim,
                               int radius,
                               int *nech_loc,
                               int *tabind,
                               double *tabval)
{
  double value;
  int i, ix, iy, iz, jpos, nech, iwork1[3], iwork2[3], nrx, nry, nrz, nmx, nmy,
      nmz;

  /* Initializations */

  nech = 0;
  nrx = (ndim >= 1) ? radius :
                      0;
  nry = (ndim >= 2) ? radius :
                      0;
  nrz = (ndim >= 3) ? radius :
                      0;
  nmx = (ndim >= 1) ? DB_GRID_FILL->getNX(0) : 1;
  nmy = (ndim >= 2) ? DB_GRID_FILL->getNX(1) : 1;
  nmz = (ndim >= 3) ? DB_GRID_FILL->getNX(2) : 1;

  /* Locate the central cell */

  for (i = 0; i < 3; i++)
    iwork1[i] = iwork2[i] = 0;
  db_index_sample_to_grid(DB_GRID_FILL, ipos, iwork1);

  /* Loop on the neighborhood cells */

  for (ix = -nrx; ix <= nrx; ix++)
  {
    iwork2[0] = iwork1[0] + ix;
    if (iwork2[0] < 0 || iwork2[0] >= nmx) continue;
    for (iy = -nry; iy <= nry; iy++)
    {
      iwork2[1] = iwork1[1] + iy;
      if (iwork2[1] < 0 || iwork2[1] >= nmy) continue;
      for (iz = -nrz; iz <= nrz; iz++)
      {
        iwork2[2] = iwork1[2] + iz;
        if (iwork2[2] < 0 || iwork2[2] >= nmz) continue;
        jpos = db_index_grid_to_sample(DB_GRID_FILL, iwork2);
        value = DB_GRID_FILL->getLocVariable(ELoc::Z,jpos, 0);
        if (FFFF(value)) continue;
        tabind[nech] = jpos;
        tabval[nech] = value;
        nech++;
      }
    }
  }

  /* Returning arguments */
  *nech_loc = nech;
  return;
}

/****************************************************************************/
/*!
 **  Calculate the extrapolation
 **
 ** \return  Error return code
 **
 ** \param[in]  ipos    Absolute grid index of the input grid node
 ** \param[in]  mode    Type of interpolation
 ** \li                 0 : Moving average
 ** \li                 1 : Inverse squared distance
 ** \li                 2 : Interpolation by a linear plane
 ** \li                 3 : Distance to the initial grains
 ** \param[in]  nech    Number of samples in the neighborhood
 ** \param[in]  tabind  Index of the neighboring sample
 ** \param[in]  tabval  Value of the neighboring sample
 **
 *****************************************************************************/
static int st_grid_fill_calculate(int ipos,
                                  int mode,
                                  int nech,
                                  int *tabind,
                                  double *tabval)
{
  double result, dist, dist2, top, bot, a[6], b[3], sol[3], f[4], dmin;
  int iech, indg[3], j, k, l, pivot;
  static int neq = 3;

  /* Initializations */

  result = top = bot = 0.;

  /* Dispatch according to the extrapolation mode */

  switch (mode)
  {
    case 0:
      for (iech = 0; iech < nech; iech++)
        top += tabval[iech];
      result = top / (double) nech;
      break;

    case 1:
      for (iech = 0; iech < nech; iech++)
      {
        dist = distance_intra(DB_GRID_FILL, ipos, tabind[iech], NULL);
        dist2 = dist * dist;
        top += tabval[iech] / dist2;
        bot += 1. / dist2;
      }
      result = top / bot;
      break;

    case 2:
      if (nech < 3) return (1);
      f[0] = 1.;
      for (iech = 0; iech < nech; iech++)
      {
        db_index_sample_to_grid(DB_GRID_FILL, tabind[iech], indg);
        grid_to_point(DB_GRID_FILL, indg, NULL, &f[1]);
        for (j = l = 0; j < neq; j++)
        {
          for (k = 0; k <= j; k++, l++)
            a[l] += f[j] * f[k];
          b[j] += tabval[iech] * f[j];
        }
      }
      if (matrix_solve(0, a, b, sol, neq, 1, &pivot)) return (1);
      if (pivot != 0) return (1);
      db_index_sample_to_grid(DB_GRID_FILL, ipos, indg);
      grid_to_point(DB_GRID_FILL, indg, NULL, &f[1]);
      for (j = 0; j < neq; j++)
        result += sol[j] * f[j];
      break;

    case 3:
      dmin = 0.;
      for (iech = 0; iech < nech; iech++)
      {
        dist = tabval[iech];
        if (dist > dmin) dmin = dist;
      }
      result = dmin + 1.;
      break;
  }

  /* Assign the result */

  DB_GRID_FILL->setLocVariable(ELoc::Z,ipos, 0, result);
  return (0);
}

/*****************************************************************************/
/*!
 **  Write a sample in the output file
 **
 ** \param[in]  db        Db characteristics
 ** \param[in]  iech      Rank of the sample
 ** \param[in]  nvar      Number of attributes
 ** \param[in]  iatt      Array of the output attribute
 ** \param[in]  tab       Array of values
 **
 *****************************************************************************/
static void st_write_active_sample(Db *db,
                                   int iech,
                                   int nvar,
                                   int *iatt,
                                   double *tab)
{
  int ivar;

  for (ivar = 0; ivar < nvar; ivar++)
    db->setArray(iech, iatt[ivar], tab[ivar]);
}

/*****************************************************************************/
/*!
 **  Check if a sample must be kept or not
 **
 ** \return  1 if the sample must be kept
 **
 ** \param[in]  db        Db characteristics
 ** \param[in]  flag_zero 1 to forbid zero values
 ** \param[in]  iech      Rank of the sample
 ** \param[in]  nvar      Number of attributes
 ** \param[in]  iatt      Array of the input attribute
 ** \param[in]  eps       Minimum value assigned to a zero
 **
 ** \param[out] tab       Array of values
 **
 *****************************************************************************/
static int st_read_active_sample(Db *db,
                                 int flag_zero,
                                 int iech,
                                 int nvar,
                                 int *iatt,
                                 double eps,
                                 double *tab)
{
  int ivar, number;

  /* Initialize the output array */

  for (ivar = 0; ivar < nvar; ivar++)
    tab[ivar] = TEST;

  /* Check for the presence of a Test value */

  for (ivar = 0; ivar < nvar; ivar++)
  {
    tab[ivar] = db->getArray(iech, iatt[ivar]);
    if (FFFF(tab[ivar])) return (0);
  }

  /* Count the number of zeroes */

  if (!flag_zero) return (1);
  number = 0;
  for (ivar = 0; ivar < nvar; ivar++)
    if (tab[ivar] <= 0.) number++;

  /* If zeroes have been found, correct values */

  if (number > 0)
  {
    for (ivar = 0; ivar < nvar; ivar++)
    {
      if (tab[ivar] <= 0.)
        tab[ivar] = eps;
      else
        tab[ivar] = tab[ivar] - number * eps / (nvar - number);
    }
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Find the interval (among vector X) to which the sample (x) belongs
 **
 ** \return  Rank of the interval containing the target
 **
 ** \param[in]  x       Target coordinate
 ** \param[in]  ndef    Number of defined samples
 ** \param[in]  X       Vector of coordinate of valued samples
 **
 ** \remarks If x < min(X) or x > max(X) the returned index is -1
 ** \remarks Array X is assumed to be increasingly ordered
 **
 *****************************************************************************/
static int st_find_interval(double x, int ndef, double *X)
{
  if (ndef < 2) return -1;
  if (x < X[0] || x > X[ndef - 1]) return -1;

  // The target belongs to an interval
  for (int k = 0; k < ndef - 1; k++)
    if (x >= X[k] && x < X[k + 1]) return k;

  // The target matches the upper bound of the last interval
  if (x == X[ndef - 1]) return ndef - 1;
  return -1;
}

/****************************************************************************/
/*!
 **  Fill an incomplete 1-D grid by linear interpolation from a set of
 **  valued samples
 **
 ** \param[in]  dbgrid  Db grid structure
 ** \param[in]  ivar    Rank of the variable to be filled
 ** \param[in]  ndef    Number of defined samples
 ** \param[in]  X       Vector of coordinate of valued samples
 ** \param[in]  Y       Vector of values of valued samples
 **
 *****************************************************************************/
static void st_grid1D_interpolate_linear(Db *dbgrid,
                                         int ivar,
                                         int ndef,
                                         double *X,
                                         double *Y)
{
  int nech = dbgrid->getSampleNumber();

  // Loop on the grid nodes

  for (int iech = 0; iech < nech; iech++)
  {
    if (!dbgrid->isActive(iech)) continue;
    double x = dbgrid->getCoordinate(iech, 0);
    int k = st_find_interval(x, ndef, X);
    if (k < 0) continue;
    double y = Y[k] + (Y[k + 1] - Y[k]) * (x - X[k]) / (X[k + 1] - X[k]);
    dbgrid->setLocVariable(ELoc::Z,iech, ivar, y);
  }
}

/****************************************************************************/
/*!
 **  Fill an incomplete 1-D grid by spline interpolation from a set of
 **  valued samples
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbgrid  Db grid structure
 ** \param[in]  ivar    Rank of the variable to be filled
 ** \param[in]  ndef    Number of defined samples
 ** \param[in]  X       Vector of coordinate of valued samples
 ** \param[in]  Y       Vector of values of valued samples
 **
 *****************************************************************************/
static int st_grid1D_interpolate_spline(Db *dbgrid,
                                        int ivar,
                                        int ndef,
                                        double *X,
                                        double *Y)
{
  VectorDouble h, F, R, M, C, Cp;
  int nech = dbgrid->getSampleNumber();

  // Preliminary calculations

  int n = ndef;
  int nm1 = n - 1;

  h.resize(nm1);
  for (int i = 0; i < nm1; i++)
    h[i] = X[i + 1] - X[i];

  F.resize(n);
  for (int i = 1; i < nm1; i++)
    F[i] = (Y[i + 1] - Y[i]) / h[i] - (Y[i] - Y[i - 1]) / h[i - 1];
  F[0] = 0;
  F[nm1] = 0;

  R.resize(n * n, 0);
  R(0,0) = 1;
  R(nm1,nm1) = 1.;
  for (int i = 1; i < nm1; i++)
  {
    R(i,i) = (h[i - 1] + h[i]) / 3;
    R(i,i+1) = h[i] / 6;
    R(i,i-1) = h[i - 1] / 6;
  }

  M.resize(n, 0);
  if (matrix_invert(R.data(), n, -1)) return 1;
  matrix_product_safe(n, n, 1, R.data(), F.data(), M.data());

  C.resize(nm1, 0);
  Cp.resize(nm1, 0);
  for (int i = 0; i < nm1; i++)
  {
    C[i] = (Y[i + 1] - Y[i]) / h[i] - (M[i + 1] - M[i]) * h[i] / 6;
    Cp[i] = Y[i] - M[i] * h[i] * h[i] / 6;
  }

  // Loop on the grid nodes

  for (int iech = 0; iech < nech; iech++)
  {
    double y = TEST;
    if (dbgrid->isActive(iech))
    {
      double x = dbgrid->getCoordinate(iech, 0);
      int k = st_find_interval(x, ndef, X);
      if (k >= 0)
      {
        double d1 = X[k + 1] - x;
        double d2 = x - X[k];
        double h6 = 6 * h[k];
        y = (M[k] * pow(d1, 3) + M[k + 1] * pow(d2, 3)) / h6 + C[k] * d2
            + Cp[k];
      }
    }
    dbgrid->setLocVariable(ELoc::Z,iech, ivar, y);
  }
  return 0;
}

/*****************************************************************************/
/*!
 **  Centers the samples of a Db to the center of blocks of a grid Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db_point   descriptor of the point parameters
 ** \param[in]  db_grid    descriptor of the grid parameters
 ** \param[in]  eps_random Randomisation Epsilon
 **
 ** \remark The argument 'eps_random' allows perturbating the centered
 ** \remark coordinate so that it does not lie exactly on the node.
 ** \remark This possibility makes sense in order to identify centered data
 ** \remark from data actually located on the grid center (before migration)
 ** \remark The perturbation is calculated as DX(i) * eps
 **
 *****************************************************************************/
int DbHelper::centerPointToGrid(Db *db_point, DbGrid *db_grid, double eps_random)
{
  if (db_point == nullptr) return 1;
  if (db_grid == nullptr) return 1;
  if (!db_point->hasSameDimension(db_grid))
  {
    messerr("For centering, 'dbin' and 'dbout' should share the same Space Dimension");
    return 1;
  }
  int ndim = db_point->getNDim();

  /* Core allocation */

  VectorDouble coor(ndim);

  /* Loop on the samples of the Point Db */

  for (int iech = 0; iech < db_point->getSampleNumber(); iech++)
  {

    /* Read the coordinates of the point sample */

    for (int idim = 0; idim < ndim; idim++)
      coor[idim] = db_point->getCoordinate(iech, idim);

    /* Get the indices of the grid node */

    db_grid->centerCoordinateInPlace(coor, true);

    /* Randomize the processed center */

    if (eps_random > 0)
      for (int idim = 0; idim < ndim; idim++)
        coor[idim] += db_grid->getDX(idim) * law_uniform(0., eps_random);

    /* Correct the sample locations */

    for (int idim = 0; idim < ndim; idim++)
      db_point->setCoordinate(iech, idim, coor[idim]);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Look for duplicates
 **
 ** \return  Error return code
 **
 ** \param[in]  db1        First Db
 ** \param[in]  db2        Second Db
 ** \param[in]  flag_same  True if the two Db files are the same
 ** \param[in]  verbose    True for verbose output
 ** \param[in]  opt_code   code selection option (if code is defined)
 ** \li                     0 : no use of the code selection
 ** \li                     1 : codes must be close enough
 ** \li                     2 : codes must be different
 ** \param[in]  tolcode    Code tolerance
 ** \param[in]  dist       Array of the minimum distance whose length is the space dimension (or NULL for a  null distance)
 **
 ** \param[out]  sel       Array containing the selection
 **
 *****************************************************************************/
int DbHelper::findDuplicates(Db *db1,
                             Db *db2,
                             bool flag_same,
                             bool verbose,
                             int opt_code,
                             double tolcode,
                             const VectorDouble &dist,
                             VectorDouble &sel)
{
  bool flag_code = db1->hasLocVariable(ELoc::C) && db2->hasLocVariable(ELoc::C);
  int nmerge = 0;

  // Title (optional)

  if (verbose) mestitle(1, "Look for duplicates");

  /* Set the selection */

  for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
    sel[iech2] = db2->getSelection(iech2);

  /* Loop on the samples of the second Db */

  for (int iech2 = 0; iech2 < db2->getSampleNumber(); iech2++)
  {
    if (!db2->isActive(iech2)) continue;
    for (int iech1 = 0; iech1 < db1->getSampleNumber(); iech1++)
    {
      if (!db1->isActive(iech1)) continue;
      if (flag_same)
      {
        if (iech1 == iech2) continue;
        if (!sel[iech1]) continue;
      }

      /* Check if the two points have similar coordinates */

      bool flag_diff = false;
      for (int idim = 0; idim < db1->getNDim() && !flag_diff; idim++)
      {
        double v1 = db1->getCoordinate(iech1, idim);
        double v2 = db2->getCoordinate(iech2, idim);
        if (flag_code)
        {
          if (st_code_comparable(db1, db2, iech1, iech2, opt_code, (int) tolcode))
            continue;
        }
        double dval = (!dist.empty()) ? dist[idim] : 0.;
        if (ABS(v1 - v2) > dval) flag_diff = true;
      }
      if (flag_diff) continue;

      sel[iech2] = 0;
      nmerge++;

      /* Optional printout */

      if (verbose)
      {
        message("Sample %d too close to sample %d\n", iech1 + 1, iech2 + 1);
        db_sample_print(db1, iech1, 1, 0, 0);
        db_sample_print(db2, iech2, 1, 0, 0);
        message("\n");
      }
    }
  }

  // Final printout (optional)

  if (verbose)
  {
    if (nmerge > 0)
      message("- Count of masked samples = %d\n", nmerge);
    else
      message("- No duplicate found\n");
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Normalize a set of variables
 **
 ** \return  Error returned code
 **
 ** \param[in]  db     Db structure
 ** \param[in]  oper   Name of the operator
 ** \li                "mean"  : Normalize the mean to the input value
 ** \li                "stdv"  : Normalize the st. dev. to the input value
 ** \li                "scal"  : Normalize the mean and st. dev.
 ** \li                "prop"  : Normalize the variables to proportions
 ** \param[in]  cols   Ranks of the variables
 ** \param[in]  center Theoretical Mean value
 ** \param[in]  stdv   Theoretical Standard Deviation value
 **
 *****************************************************************************/
int DbHelper::normalizeVariables(Db *db,
                                 const char *oper,
                                 const VectorInt& cols,
                                 double center,
                                 double stdv)
{
  int jcol, ndef, iptr;
  double *num, *mm, *vv, proptot, value;

  /* Initializations */

  int nech = db->getSampleNumber();
  num = mm = vv = nullptr;
  int ncol = (int) cols.size();

  /* Check that all variables are defined */

  for (int icol = 0; icol < ncol; icol++)
  {
    jcol = cols[icol];
    if (!db->isColIdxValid(jcol))
    {
      messerr("Column %d is not defined", cols[icol]);
      return (1);
    }
  }

  /* Core allocation */

  num = (double*) mem_alloc(sizeof(double) * ncol, 0);
  if (num == nullptr) goto label_end;
  mm = (double*) mem_alloc(sizeof(double) * ncol, 0);
  if (mm == nullptr) goto label_end;
  vv = (double*) mem_alloc(sizeof(double) * ncol, 0);
  if (vv == nullptr) goto label_end;

  /* Initializations */

  for (int icol = 0; icol < ncol; icol++)
    num[icol] = mm[icol] = vv[icol] = 0.;

  /* Printout */

  if (!strcmp(oper, "mean"))
  {
    if (FFFF(center)) center = 0.;
  }
  else if (!strcmp(oper, "stdv"))
  {
    if (FFFF(stdv)) stdv = 1.;
  }
  else if (!strcmp(oper, "scal"))
  {
    if (FFFF(center)) center = 0.;
    if (FFFF(stdv)) stdv = 1.;
  }
  else if (!strcmp(oper, "prop"))
  {
    center = 0.;                  // Useless line
  }
  else
  {
    messerr("Invalid operator name (%s)", oper);
    messerr("The list of operators available is:");
    messerr("mean  : Normalize the mean");
    messerr("stdv  : Normalize the st. dev.");
    messerr("scal  : Normalize the mean and st dev");
    messerr("prop  : Normalize the proportions");
    return (1);
  }

  /* Creation of the new attributes */

  iptr = db->addColumnsByConstant(ncol, TEST);
  if (iptr < 0) return (1);

  /* Loop on the samples */

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Loop on the variables */

    ndef = 0;
    proptot = 0.;
    for (int icol = 0; icol < ncol; icol++)
    {
      value = db->getArray(iech, cols[icol]);
      if (FFFF(value)) continue;
      if (!strcmp(oper, "prop")) value = MIN(1., MAX(0.,value));

      /* Update statistics */

      ndef += 1;
      num[icol] += 1;
      proptot += value;
      mm[icol] += value;
      vv[icol] += value * value;
    }

    /* Set the output array (in the case of proportions) */

    if (!strcmp(oper, "prop"))
    {
      for (int icol = 0; icol < ncol; icol++)
      {
        value = db->getArray(iech, cols[icol]);
        value = MIN(1., MAX(0.,value));
        if (ndef == ncol && proptot > 0.)
          db->setArray(iech, iptr + icol, value / proptot);
      }
    }
  }

  if (strcmp(oper, "prop"))
  {

    /* Global Normation */

    for (int icol = 0; icol < ncol; icol++)
    {
      if (num[icol] <= 0)
      {
        mm[icol] = TEST;
        vv[icol] = TEST;
      }
      else
      {
        mm[icol] /= num[icol];
        vv[icol] = vv[icol] / num[icol] - mm[icol] * mm[icol];
      }
    }

    /* Set the output array */

    for (int iech = 0; iech < nech; iech++)
    {
      for (int icol = 0; icol < ncol; icol++)
      {
        jcol = cols[icol];
        value = db->getArray(iech, jcol);
        if (!strcmp(oper, "mean"))
        {
          if (!FFFF(mm[icol]))
            db->setArray(iech, iptr + icol, value + center - mm[icol]);
        }
        else if (!strcmp(oper, "stdv"))
        {
          if (!FFFF(vv[icol]) && vv[icol] > 0.)
            db->setArray(iech, iptr + icol, value * stdv / sqrt(vv[icol]));
        }
        else if (!strcmp(oper, "scal"))
        {
          if (!FFFF(vv[icol]) && vv[icol] > 0. && !FFFF(mm[icol]))
            db->setArray(iech, iptr + icol,
                         center + (value - mm[icol]) * stdv / sqrt(vv[icol]));
        }
      }
    }
  }

  label_end: num = (double*) mem_free((char* ) num);
  mm = (double*) mem_free((char* ) mm);
  vv = (double*) mem_free((char* ) vv);

  return (0);
}

/****************************************************************************/
/*!
 **  Fill an incomplete grid
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbgrid  Db grid structure
 ** \param[in]  mode    Type of interpolation
 ** \li                  0 : Moving average
 ** \li                  1 : Inverse squared distance
 ** \li                  2 : Interpolation by a linear plane
 ** \li                  3 : Distance to the initial grains
 ** \param[in]  seed    Seed used for the random number generation
 ** \param[in]  radius  Radius of the neighborhood
 ** \param[in]  verbose Verbose flag
 ** \param[in]  namconv Naming convention
 **
 *****************************************************************************/
int DbHelper::dbgrid_filling(DbGrid *dbgrid,
                             int mode,
                             int seed,
                             int radius,
                             bool verbose,
                             const NamingConvention &namconv)
{
  Skin *skin = nullptr;
  double *tabval;
  int *tabind, error, rank, ipos, ndim, count, nech;
  LocalSkin SKF;

  /* Initializations */

  error = 1;

  /* Preliminary checks */

  if (! dbgrid->isGrid())
  {
    messerr("This function is limited to Grid Db");
    return (1);
  }
  if (!dbgrid->isVariableNumberComparedTo(1)) return (1);
  ndim = dbgrid->getNDim();
  if (ndim > 3)
  {
    messerr("This function is limited to a maximum 3-D space");
    return (1);
  }
  if (mode < 0 || mode > 3)
  {
    messerr("The argument 'mode' should lie between 0 and 3");
    return (1);
  }
  if (mode == 2 && radius < 1)
  {
    messerr("The linear interpolation requires a neighborhood radius > 1");
    return (1);
  }

  // Create the new variable and duplicate the Z-locator variable

  int iatt_in = dbgrid->getUIDByLocator(ELoc::Z, 0);
  int iatt_out = dbgrid->addColumnsByConstant(1);
  dbgrid->duplicateColumnByUID(iatt_in, iatt_out);
  dbgrid->setLocatorByUID(iatt_out, ELoc::Z);

  /* Global variables */

  DB_GRID_FILL = dbgrid;
  skin = nullptr;
  tabval = nullptr;
  tabind = nullptr;
  count = (int) pow(2. * radius + 1., (double) ndim) - 1;

  /* Core allocation */

  law_set_random_seed(seed);
  tabval = (double*) mem_alloc(sizeof(double) * count, 0);
  if (tabval == nullptr) goto label_end;
  tabind = (int*) mem_alloc(sizeof(int) * count, 0);
  if (tabind == nullptr) goto label_end;


  skin = new Skin(&SKF, dbgrid);

  if (skin->init(verbose))
  {
    error = 0;
    goto label_end;
  }

  /* Implicit loop on the cells to be filled */

  while (skin->remains(verbose))
  {

    /* Find the next cell to be processed */

    skin->getNext(&rank, &ipos);

    /* Find the neighborhood */

    st_grid_fill_neigh(ipos, ndim, radius, &nech, tabind, tabval);

    /* Calculate the extrapolated value */

    if (st_grid_fill_calculate(ipos, mode, nech, tabind, tabval)) continue;

    /* Deduce the initial influence of the central cell */

    if (skin->unstack(rank, ipos)) goto label_end;
  }

  // Optional printout

  if (verbose) skin->skinPrint();

  /* Set the error return code */

  error = 0;
  namconv.setNamesAndLocators(dbgrid, {iatt_in}, dbgrid, iatt_out);

  label_end:
  delete skin;
  tabval = (double*) mem_free((char* ) tabval);
  tabind = (int*) mem_free((char* ) tabind);
  return (error);
}

/****************************************************************************/
/*!
 **  Look for duplicates within a Db
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db Structure
 ** \param[in]  verbose    True for verbose output
 ** \param[in]  dist       Array of the minimum distance
 ** \param[in]  opt_code   code selection option (if code is defined)
 ** \li                     0 : no use of the code selection
 ** \li                     1 : codes must be close enough
 ** \li                     2 : codes must be different
 ** \param[in]  tolcode    Code tolerance
 ** \param[in]  namconv    Naming convention
 **
 *****************************************************************************/
int DbHelper::db_duplicate(Db *db,
                           bool verbose,
                           const VectorDouble &dist,
                           int opt_code,
                           double tolcode,
                           const NamingConvention &namconv)
{
  if (db == nullptr)
  {
    messerr("You must define a Db");
    return 1;
  }

  // Adding a new variable

  VectorDouble sel(db->getSampleNumber());

  // Check for duplicates

  if (DbHelper::findDuplicates(db, db, 1, verbose, opt_code, tolcode, dist,
                               sel)) return 1;

  // Add the variable to the Db
  int iatt = db->addColumns(sel);

  // Setting the output variable
  namconv.setNamesAndLocators(db, iatt);

  return 0;
}

/*****************************************************************************/
/*!
 **  Translate a set of compositional variables into auxiliary variables
 **
 ** \return  Error return code
 **
 ** \param[in]  db       Db characteristics
 ** \param[in]  verbose  1 for a Verbose option
 ** \param[in]  mode     1 for Forward and -1 for Backward transformation
 ** \param[in]  type     Type of conversion
 ** \li                  0 : Simple transformation in proportion
 ** \li                  1 : Additive logratio
 ** \li                  2 : Centered logratio
 ** \li                  3 : Isometric logratio
 ** \param[in]  number   Number of input attributes
 ** \param[in]  iatt_in  Array of the input attribute
 ** \param[in]  iatt_out Array of the output attribute
 **
 ** \param[out] numout   Number of variables in output
 **
 ** \remarks  The additive and the isometric logratio transformations
 ** \remarks  transform N+1 compositional variables into N elements (Forward)
 ** \remarks  and from N transformed elements into N+1 compositional variables
 ** \remarks  (Backwards).
 **
 ** \remarks  The zero-values are replaced by a conventional small value
 ** \remarks  This is defined by a variable that can be corrected using
 ** \remarks  the keypair facility with the keyword 'CompositionalEps'
 **
 ** \remarks  The arguments 'iatt_in' and 'iatt_out' can coincide.
 **
 ** \remarks  Outlier Detection for Compositional Data Using Robust Methods
 ** \remarks  Math Geosciences (2008) 40: 233-248
 **
 *****************************************************************************/
int DbHelper::db_compositional_transform(Db *db,
                                         int verbose,
                                         int mode,
                                         int type,
                                         int number,
                                         int *iatt_in,
                                         int *iatt_out,
                                         int *numout)
{
  int error, nech, number1, iech, ivar, jvar;
  double *tabin, *tabout, sum, eps;

  /* Initializations */

  error = 1;
  nech = db->getSampleNumber();
  tabin = tabout = nullptr;
  eps = get_keypone("CompositionalEps", EPSILON3);

  /* Core allocation (may be one more than needed, but general) */

  number1 = number + 1;
  tabin = (double*) mem_alloc(sizeof(double) * number1, 0);
  if (tabin == nullptr) goto label_end;
  tabout = (double*) mem_alloc(sizeof(double) * number1, 0);
  if (tabout == nullptr) goto label_end;

  /* Verbose output */

  if (verbose) mestitle(0, "Compositional transformation");

  /* Dispatch */

  switch (type)
  {
    case 0: /* No action */
      if (verbose)
      {
        if (mode > 0)
          message("- Scaling to Proportions (Forwards)\n");
        else
          message("- Scaling to Proportions (Backwards)\n");
      }
      (*numout) = number;
      for (iech = 0; iech < nech; iech++)
      {
        if (!db->isActive(iech)) continue;
        (void) st_read_active_sample(db, 1, iech, number, iatt_in, eps, tabin);
        sum = 0.;
        for (ivar = 0; ivar < (*numout); ivar++)
          sum += tabin[ivar];
        for (ivar = 0; ivar < (*numout); ivar++)
          tabout[ivar] = (sum > 0) ? tabin[ivar] / sum :
                                     TEST;
        st_write_active_sample(db, iech, *numout, iatt_out, tabout);
      }
      break;

    case 1: /* Additive Logratio */
      if (mode > 0)
      {
        /* Forward transformation */

        if (verbose) message("- using Additive Log-Ratio (Forwards)\n");
        (*numout) = number - 1;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 1, iech, number, iatt_in, eps, tabin))
          {
            for (ivar = 0; ivar < (*numout); ivar++)
              tabout[ivar] = log(tabin[ivar] / tabin[(*numout)]);
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      else
      {
        /* Backward transformation */

        message("- using Additive Log-Ratio (Backwards)\n");
        (*numout) = number + 1;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 0, iech, number, iatt_in, eps, tabin))
          {
            sum = 1.;
            for (ivar = 0; ivar < number; ivar++)
            {
              tabout[ivar] = exp(tabin[ivar]);
              sum += tabout[ivar];
            }
            for (ivar = 0; ivar < number; ivar++)
              tabout[ivar] /= sum;
            tabout[number] = 1 / sum;
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      break;

    case 2: /* Centered Logratio */
      if (mode > 0)
      {
        /* Forward transformation */

        if (verbose) message("- using Centered Log-Ratio (Forwards)\n");
        (*numout) = number;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 1, iech, number, iatt_in, eps, tabin))
          {
            sum = 0.;
            for (ivar = 0; ivar < number; ivar++)
            {
              tabout[ivar] = log(tabin[ivar]);
              sum += tabout[ivar];
            }
            for (ivar = 0; ivar < number; ivar++)
              tabout[ivar] = tabout[ivar] - 0.5 * sum;
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      else
      {
        /* Backward transformation */

        if (verbose) message("- using Centered Log-Ratio (Backwards)\n");
        (*numout) = number;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 0, iech, number, iatt_in, eps, tabin))
          {
            sum = 0.;
            for (ivar = 0; ivar < number; ivar++)
            {
              tabout[ivar] = exp(tabin[ivar]);
              sum += tabout[ivar];
            }
            for (ivar = 0; ivar < number; ivar++)
              tabout[ivar] /= sum;
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      break;

    case 3: /* Isometric Logratio */
      if (mode > 0)
      {
        /* Forward transformation */

        if (verbose) message("- using Isometric Log-Ratio (Forwards)\n");
        (*numout) = number - 1;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 1, iech, number, iatt_in, eps, tabin))
          {
            for (ivar = 0; ivar < number; ivar++)
              tabin[ivar] = log(tabin[ivar]);
            sum = 0.;
            for (ivar = 0; ivar < (*numout); ivar++)
            {
              sum += tabin[ivar];
              tabout[ivar] = (sum / (ivar + 1.) - tabin[ivar + 1])
                  * sqrt((ivar + 1.) / (ivar + 2.));
            }
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      else
      {
        /* Backward transformation */

        if (verbose) message("- using Isometric Log-Ratio (Backwards)\n");
        (*numout) = number + 1;
        for (iech = 0; iech < nech; iech++)
        {
          if (!db->isActive(iech)) continue;
          if (st_read_active_sample(db, 0, iech, number, iatt_in, eps, tabin))
          {
            tabin[number] = 0.;
            for (ivar = 0; ivar < (*numout); ivar++)
            {
              sum = 0.;
              for (jvar = ivar; jvar < (*numout); jvar++)
                sum += tabin[jvar] / sqrt((jvar + 1.) * (jvar + 2.));
              tabout[ivar] = sum;
              if (ivar > 0)
                tabout[ivar] -= tabin[ivar - 1] * sqrt(ivar / (ivar + 1.));
            }

            sum = 0.;
            for (jvar = 0; jvar < (*numout); jvar++)
            {
              tabout[jvar] = exp(tabout[jvar]);
              sum += tabout[jvar];
            }
            for (ivar = 0; ivar < (*numout); ivar++)
              tabout[ivar] /= sum;
          }
          st_write_active_sample(db, iech, *numout, iatt_out, tabout);
        }
      }
      break;
  }

  /* Set the error return code */

  error = 0;

  label_end: tabin = (double*) mem_free((char* ) tabin);
  tabout = (double*) mem_free((char* ) tabout);
  return (error);
}

/*****************************************************************************/
/*!
 **  Sample a grid into a finer subgrid (all variables)
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin  Descriptor of the grid parameters
 ** \param[in]  nmult Array of multiplicity coefficients
 **
 *****************************************************************************/
DbGrid* DbHelper::dbgrid_sampling(DbGrid *dbin, const VectorInt &nmult)
{
  DbGrid *dbout;
  VectorDouble coor;
  int ncol, icol, iech, iad, item, rank, ndim;
  ELoc locatorType;

  /* Initializations */

  dbout = nullptr;
  ncol = dbin->getColumnNumber();
  ndim = dbin->getNDim();

  /* Core allocation */

  coor.resize(ndim);

  /* Create the subgrid */

  dbout = db_create_grid_multiple(dbin, nmult, 1);
  if (dbout == nullptr) goto label_end;
  rank = dbout->addColumnsByConstant(ncol, TEST);
  if (rank < 0) goto label_end;
  for (icol = 0; icol < ncol; icol++)
  {
    (void) dbin->getLocatorByColIdx(icol, &locatorType, &item);
    dbout->setLocatorByUID(icol, locatorType, item);
  }

  /* Loop on the samples of the output grid */

  for (iech = 0; iech < dbout->getSampleNumber(); iech++)
  {
    if (!dbout->isActive(iech)) continue;
    db_sample_load(dbout, ELoc::X, iech, coor.data());
    iad = dbin->coordinateToRank(coor);
    if (iad < 0) continue;

    /* Loop on the variables of the input grid */

    for (icol = 0; icol < ncol; icol++)
      dbout->setValueByColIdx(iech, icol, dbin->getValueByColIdx(iad, icol));
  }

  label_end:
  return (dbout);
}

/****************************************************************************/
/*!
 **  Fill an incomplete 1-D grid
 **
 ** \return  Error returned code
 **
 ** \param[in]  dbgrid  Db grid structure
 ** \param[in]  mode    Type of interpolation
 ** \li                  0 : Linear interpolation
 ** \li                  1 : Cubic Spline
 ** \param[in]  seed    Seed used for the random number generation
 ** \param[in]  namconv Naming convention
 **
 *****************************************************************************/
int DbHelper::db_grid1D_fill(DbGrid *dbgrid,
                             int mode,
                             int seed,
                             const NamingConvention &namconv)
{
  /* Preliminary checks */

  if (! dbgrid->isGrid())
  {
    messerr("This function is limited to Grid Db");
    return (1);
  }
  int ndim = dbgrid->getNDim();
  if (ndim != 1)
  {
    messerr("This function is limited to 1-D space");
    return (1);
  }
  if (mode < 0 || mode > 1)
  {
    messerr("The argument 'mode' should lie between 0 and 1");
    return (1);
  }
  int nvar = dbgrid->getLocNumber(ELoc::Z);
  if (nvar <= 0)
  {
    messerr("You must have at least one Z-locator defined");
    return 1;
  }

  // Add the variables (they must be defined as ELoc::Z) for following functions

  int iatt_out = dbgrid->addColumnsByConstant(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int iatt_in = dbgrid->getUIDByLocator(ELoc::Z, ivar);
    dbgrid->duplicateColumnByUID(iatt_in, iatt_out + ivar);
  }

  int nech = dbgrid->getSampleNumber();

  /* Core allocation */

  law_set_random_seed(seed);
  VectorDouble X(nech, 0);
  VectorDouble Y(nech, 0);

  /* Loop on the variables to be filled */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    // Copy the input variable

    // Look for the defined values

    int ndef = 0;
    for (int iech = 0; iech < nech; iech++)
    {
      if (!dbgrid->isActive(iech)) continue;
      double value = dbgrid->getLocVariable(ELoc::Z,iech, ivar);
      if (FFFF(value)) continue;
      X[ndef] = dbgrid->getCoordinate(iech, 0);
      Y[ndef] = value;
      ndef++;
    }

    // Perform the 1-D interpolation

    if (mode == 0)
      st_grid1D_interpolate_linear(dbgrid, ivar, ndef, X.data(), Y.data());
    else
    {
      if (st_grid1D_interpolate_spline(dbgrid, ivar, ndef, X.data(), Y.data()))
        return 1;
    }
  }

  /* Set the error return code */

  namconv.setNamesAndLocators(dbgrid, VectorString(), ELoc::Z, -1, dbgrid, iatt_out);

  return 0;
}

