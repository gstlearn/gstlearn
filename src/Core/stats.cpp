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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Morpho/Morpho.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Basic/EJustify.hpp"
#include "Basic/File.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"

#include <math.h>
#include <string.h>

/*! \cond */
#define G_ADDRESS(ix,iy,iz,nxyz)    ((ix) + nxyz[0] * ((iy) + nxyz[1] * (iz)))
#define N1_TAB(ix,iy,iz) (numtab1[G_ADDRESS(ix,iy,iz,nxyz1)])
#define N2_TAB(ix,iy,iz) (numtab2[G_ADDRESS(ix,iy,iz,nxyz2)])
#define V1_TAB(ix,iy,iz) (valtab1[G_ADDRESS(ix,iy,iz,nxyz1)])
#define V2_TAB(ix,iy,iz) (valtab2[G_ADDRESS(ix,iy,iz,nxyz2)])
#define D1_TAB(ix,iy,iz) (N1_TAB(ix,iy,iz) > 0 &&     \
                          ! FFFF(V1_TAB(ix,iy,iz)) && \
                          V1_TAB(ix,iy,iz) > 0)
#define RESIDUALS(icut,iech) (residuals[(icut) * nech + (iech)])
#define NBGH(ivois,idim)     (nbgh[ndim * (ivois) + (idim)])
#define TABINI(iseed,idim)   (tabini[ndim * (iseed) + (idim)])
#define TABCUR(iseed,idim)   (tabcur[ndim * (iseed) + (idim)])
#define TRAJEC(iseed,iter,idim) (trsave[(niter * (iseed) + (iter)) * ndim + (idim)])
/*! \endcond */

static int DEBUG = 0;

/****************************************************************************/
/*!
 **  Copy the multivariate into monovariate statistics (before printout)
 **
 ** \param[in]  ncol  Dimension of the (square) matrix
 ** \param[in,out]  tab   Array to be refactored
 **
 *****************************************************************************/
static void st_refactor(int ncol, double *tab)
{
  int ix, iy;

  for (ix = 0; ix < ncol; ix++)
    for (iy = 0; iy < ncol; iy++)
      tab[ix * ncol + iy] = tab[iy * ncol + iy];
}

/****************************************************************************/
/*!
 **  Copy the multivariate or monovariate statistics into the returned array
 **
 ** \param[in]  nx    First dimension of the matrix
 ** \param[in]  ny    First dimension of the matrix
 ** \param[in]  tabin Array to be refactored
 **
 ** \param[out]  tabout Array to be refactored
 **
 *****************************************************************************/
static void st_copy_result(int nx, int ny, double *tabin, double *tabout)
{
  int ix, iy, lec;

  for (ix = lec = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++, lec++)
      tabout[lec] = tabin[lec];
}

/****************************************************************************/
/*!
 **  Check the operator name is mentionned within a list
 **
 ** \return  1 if the operator is mentionned; 0 otherwise
 **
 ** \param[in]  opers Array of operators
 ** \param[in]  refe  Reference operator
 **
 ** \remarks If the array 'opers' if empty, any name is considered as valid
 **
 *****************************************************************************/
static int st_oper_exists(const VectorString &opers, const String &refe)
{
  int noper = static_cast<int>(opers.size());
  if (noper == 0) return (1);
  for (int i = 0; i < noper; i++)
  {
    if (opers[i] == refe) return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Check the operator name
 **
 ** \return  1 if the operator is valid; 0 otherwise
 **
 ** \param[in]  oper Name of the operator
 ** \param[in]  flag_multi  1 if multivariate operator is authorized
 ** \param[in]  flag_indic  1 if indicator ("plus","minus","zero") is authorized
 ** \param[in]  flag_sum    1 if sum of variable is authorized
 ** \param[in]  flag_median 1 if median is quathorized
 ** \param[in]  flag_qt     1 if QT ("ore","metal") is authorized
 **
 ** \remarks If an error occurred, the message is printed
 **
 *****************************************************************************/
static int st_oper_check(const String &oper,
                         int flag_multi,
                         int flag_indic,
                         int flag_sum,
                         int flag_median,
                         int flag_qt)
{
  int valid;

  /* Initializations */

  valid = 0;

  /* Monovariate check */

  if (oper == "num") valid = 1;
  if (oper == "mean") valid = 1;
  if (oper == "var") valid = 1;
  if (oper == "corr") valid = 1;
  if (oper == "stdv") valid = 1;
  if (oper == "mini") valid = 1;
  if (oper == "maxi") valid = 1;
  if (flag_sum)
  {
    if (oper == "sum") valid = 1;
  }
  if (flag_median)
  {
    if (oper == "med") valid = 1;
  }

  /* Multivariate check */

  if (flag_multi)
  {
    if (oper == "mean2") valid = 1;
    if (oper == "var2") valid = 1;
    if (oper == "stdv2") valid = 1;
    if (flag_sum)
    {
      if (oper == "sum2") valid = 1;
    }
  }

  /* Indicator check */

  if (flag_indic)
  {
    if (oper == "plus") valid = 1;
    if (oper == "moins") valid = 1;
    if (oper == "zero") valid = 1;
  }

  /* QT variables check */

  if (flag_qt)
  {
    if (oper == "ore") valid = 1;
    if (oper == "metal") valid = 1;
  }

  if (!valid)
  {
    messerr("Invalid operator name: '%s'", oper.c_str());
    messerr("The available operators are:");

    /* Monovariate operators */

    messerr("num   : Number of defined values");
    messerr("mean  : Mean over the defined values");
    messerr("var   : Variance over the defined values");
    messerr("corr  : Correlation over the defined values");
    messerr("stdv  : Standard Deviation over the defined values");
    messerr("mini  : Minimum over the defined values");
    messerr("maxi  : Maximum over the defined values");
    if (flag_sum)
    {
      messerr("sum   : Sum of the defined values");
    }

    /* Multivariate operators */

    if (flag_multi)
    {
      messerr("mean2 : Mean of the second variable");
      messerr("sum2  : Sum of the second variable");
      messerr("stdv2 : Standard deviation of the second variable");
      messerr("var2  : Variance of the second variable");
      if (flag_sum)
      {
        messerr("sum2  : Sum of the second variable");
      }
    }

    /* Indicator variables */

    if (flag_indic)
    {
      messerr("plus  : Number of positive samples");
      messerr("moins : Number of negative samples");
      messerr("zero  : Number of zero samples");
    }

    /* Qt variables */

    if (flag_qt)
    {
      messerr("ore   : Tonnage quantity above each cutoff");
      messerr("metal : Metal quantity above each cutoff");
    }
  }

  return (valid);
}

/****************************************************************************/
/*!
 **  Calculate the multivariate statistics between different variables of a Db
 **
 ** \return  Error Return code
 **
 ** \param[in]  db   Db structure
 ** \param[in]  oper Name of the operator
 ** \li              "num"   : Number of valid samples
 ** \li              "mean"  : Mean over the valid samples
 ** \li              "var"   : Variance over the valid samples
 ** \li              "corr"  : Correlation over the valid samples
 ** \li              "stdv"  : Standard Deviation over the valid samples
 ** \li              "mini"  : Minimum over the valid samples
 ** \li              "maxi"  : Maximum over the valid samples
 ** \li              "plus"  : Number of positive values
 ** \li              "moins" : Number of negative values
 ** \li              "zero"  : Number of zero values
 ** \param[in]  cols Ranks of the variables
 ** \param[in]  flag_mono 1 for a monovariate output
 ** \param[in]  flag_verbose  1 for a verbose output
 **
 ** \param[out] result    Resulting array
 **                       (Dimension: ncol (if flag_mono) or ncol*ncol)
 **
 *****************************************************************************/
int db_stats(Db *db,
                             const String &oper,
                             const VectorInt &cols,
                             int flag_mono,
                             int flag_verbose,
                             double *result)
{
  int i, iech, nech, icol, jcol, icol1, jcol1, icol2, jcol2, iad, ncol2;
  int error, nx, ny, ncol;
  double *num, *m1, *m2, *v1, *v2, *v12, *mini, *maxi, *plus, *moins, *zero;
  double val1, val2, weight;

  /* Initializations */
  error = 1;
  nech = db->getSampleNumber();
  ncol = static_cast<int>(cols.size());
  ncol2 = ncol * ncol;
  num = m1 = m2 = v1 = v2 = v12 = mini = maxi = nullptr;
  plus = moins = zero = nullptr;

  /* Check that all variables are defined */

  for (icol = 0; icol < ncol; icol++)
  {
    jcol = cols[icol];
    if (!db->isColumnIndexValid(jcol))
    {
      messerr("Error: Variable %d is not defined", cols[icol]);
      goto label_end;
    }
  }

  /* Check the validity of the operator */

  if (!st_oper_check(oper, 0, 1, 0, 0, 0)) goto label_end;

  /* Core allocation */

  num = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (num == nullptr) goto label_end;
  m1 = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (m1 == nullptr) goto label_end;
  m2 = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (m2 == nullptr) goto label_end;
  v1 = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (v1 == nullptr) goto label_end;
  v2 = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (v2 == nullptr) goto label_end;
  v12 = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (v12 == nullptr) goto label_end;
  mini = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (mini == nullptr) goto label_end;
  maxi = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (maxi == nullptr) goto label_end;
  plus = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (plus == nullptr) goto label_end;
  moins = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (moins == nullptr) goto label_end;
  zero = (double*) mem_alloc(sizeof(double) * ncol2, 0);
  if (zero == nullptr) goto label_end;

  /* Initializations */

  for (i = 0; i < ncol * ncol; i++)
  {
    num[i] = m1[i] = m2[i] = v1[i] = v2[i] = v12[i] = 0.;
    plus[i] = moins[i] = zero[i] = 0;
    mini[i] = 1.e30;
    maxi[i] = -1.e30;
  }

  /* Loop on the samples */

  for (iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;
    weight = db->getWeight(iech);

    /* Loop on the first variable */

    for (icol1 = 0; icol1 < ncol; icol1++)
    {
      jcol1 = cols[icol1];
      val1 = db->getArray(iech, jcol1);
      if (FFFF(val1)) continue;

      /* Loop on the second variable */

      for (icol2 = 0; icol2 < ncol; icol2++)
      {
        jcol2 = cols[icol2];
        val2 = db->getArray(iech, jcol2);
        if (FFFF(val2)) continue;

        /* Update statistics */

        iad = icol1 * ncol + icol2;
        num[iad] += weight;
        m1[iad] += weight * val1;
        m2[iad] += weight * val2;
        v1[iad] += weight * val1 * val1;
        v2[iad] += weight * val2 * val2;
        v12[iad] += weight * val1 * val2;
        if (val1 < mini[iad]) mini[iad] = val1;
        if (val1 > maxi[iad]) maxi[iad] = val1;
        if (val1 > 0) plus[iad] += 1;
        if (val1 < 0) moins[iad] += 1;
        if (val1 == 0) zero[iad] += 1;
      }
    }
  }

  /* Normation */

  for (icol1 = 0; icol1 < ncol; icol1++)
    for (icol2 = 0; icol2 < ncol; icol2++)
    {
      iad = icol1 * ncol + icol2;
      if (num[iad] <= 0)
      {
        m1[iad] = TEST;
        m2[iad] = TEST;
        v1[iad] = TEST;
        m2[iad] = TEST;
        v12[iad] = TEST;
        mini[iad] = TEST;
        maxi[iad] = TEST;
        plus[iad] = TEST;
        moins[iad] = TEST;
        zero[iad] = TEST;
      }
      else
      {
        m1[iad] /= num[iad];
        m2[iad] /= num[iad];
        v1[iad] = v1[iad] / num[iad] - m1[iad] * m1[iad];
        v2[iad] = v2[iad] / num[iad] - m2[iad] * m2[iad];
        v12[iad] = v12[iad] / num[iad] - m1[iad] * m2[iad];

        if (oper == "stdv")
        {
          v12[iad] = (v12[iad] > 0) ? sqrt(v12[iad]) :
                                      0.;
        }
        if (oper == "corr")
        {
          v12[iad] =
              (v1[iad] > 0 && v2[iad] > 0) ? v12[iad] / sqrt(
                                                 v1[iad] * v2[iad]) :
                                             0.;
        }
      }
    }

  /* Printout */

  if (flag_mono)
  {
    nx = 1;
    ny = ncol;
    st_refactor(ncol, num);
    st_refactor(ncol, m1);
    st_refactor(ncol, v1);
    st_refactor(ncol, v12);
    st_refactor(ncol, mini);
    st_refactor(ncol, maxi);
    st_refactor(ncol, plus);
    st_refactor(ncol, moins);
    st_refactor(ncol, zero);
  }
  else
  {
    nx = ncol;
    ny = ncol;
  }

  /* Optional printout */

  if (flag_verbose)
  {
    message("\n");
    if (oper == "num")
      print_matrix("Matrix of Number of defined samples", 0, 1, nx, ny, NULL,
                   num);
    else if (oper == "mean")
      print_matrix("Matrix of Variable Means", 0, 1, nx, ny, NULL, m1);
    else if (oper == "var")
      print_matrix("Matrix of Variable Variances", 0, 1, nx, ny, NULL, v12);
    else if (oper == "corr")
      print_matrix("Matrix of Variable Correlations", 0, 1, nx, ny, NULL, v12);
    else if (oper == "stdv")
      print_matrix("Matrix of Variable Standard Deviations", 0, 1, nx, ny, NULL,
                   v12);
    else if (oper == "mini")
      print_matrix("Matrix of Variable Minima", 0, 1, nx, ny, NULL, mini);
    else if (oper == "maxi")
      print_matrix("Matrix of Variable Maxima", 0, 1, nx, ny, NULL, maxi);
    else if (oper == "plus")
      print_matrix("Matrix of Number of positive samples", 0, 1, nx, ny, NULL,
                   plus);
    else if (oper == "moins")
      print_matrix("Matrix of Number of negative samples", 0, 1, nx, ny, NULL,
                   moins);
    else if (oper == "zero")
      print_matrix("Matrix of Number of zero samples", 0, 1, nx, ny, NULL,
                   zero);
    else
      messageAbort("This error should never happen");
  }

  /* Set the return array */

  if (oper == "num")
    st_copy_result(nx, ny, num, result);
  else if (oper == "mean")
    st_copy_result(nx, ny, m1, result);
  else if (oper == "var")
    st_copy_result(nx, ny, v12, result);
  else if (oper == "corr")
    st_copy_result(nx, ny, v12, result);
  else if (oper == "stdv")
    st_copy_result(nx, ny, v12, result);
  else if (oper == "mini")
    st_copy_result(nx, ny, mini, result);
  else if (oper == "maxi")
    st_copy_result(nx, ny, maxi, result);
  else if (oper == "plus")
    st_copy_result(nx, ny, plus, result);
  else if (oper == "moins")
    st_copy_result(nx, ny, moins, result);
  else if (oper == "zero")
    st_copy_result(nx, ny, zero, result);
  else
    messageAbort("This error should never happen");

  /* Set the error return code */

  error = 0;

  label_end: num = (double*) mem_free((char* ) num);
  m1 = (double*) mem_free((char* ) m1);
  m2 = (double*) mem_free((char* ) m2);
  v1 = (double*) mem_free((char* ) v1);
  v2 = (double*) mem_free((char* ) v2);
  v12 = (double*) mem_free((char* ) v12);

  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the indices of the cell neighboring a target cell
 **
 ** \param[in]  ndim   Space dimension
 ** \param[in]  radius Neighborhood radius
 ** \param[in]  rank0  Rank of the neighbor
 ** \param[in]  indg0  Array of indices of the target cell
 **
 ** \param[out] indg   Array of indices of the neighboring cell
 **
 *****************************************************************************/
static void st_get_neighboring_cell(int ndim,
                                    int radius,
                                    int rank0,
                                    int *indg0,
                                    int *indg)
{
  int nei1d, value, divid, count, reste, ratio, idim;

  /* Initializtions */

  nei1d = 2 * radius + 1;
  count = (int) pow(nei1d, (double) ndim);
  value = rank0;

  /* Loop on the space dimensions */

  divid = count;
  for (int jdim = 0; jdim < ndim; jdim++)
  {
    idim = ndim - jdim - 1;
    divid /= nei1d;
    ratio = value / divid;
    reste = value - ratio * divid;
    value = reste;
    indg[idim] = indg0[idim] + ratio - radius;
  }
}

/****************************************************************************/
/*!
 **  Calculates the monovariate statistics within cells of a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  db     Db for the points
 ** \param[in]  dbgrid Db for the grid
 ** \param[in]  oper   Name of the operator
 ** \li                "num"   : Number of valid samples
 ** \li                "mean"  : Mean over the valid samples
 ** \li                "var"   : Variance over the valid samples
 ** \li                "stdv"  : Standard Deviation over the valid samples
 ** \li                "mini"  : Minimum over the valid samples
 ** \li                "maxi"  : Maximum over the valid samples
 ** \li                "med"   : Median over the valid samples
 ** \li                "plus"  : Number of positive values
 ** \li                "moins" : Number of negative values
 ** \li                "zero"  : Number of zero values
 ** \param[in]  ncol   Number of variables
 ** \param[in]  cols   Ranks of the variables
 ** \param[in]  radius Neighborhood radius
 **
 *****************************************************************************/
int db_stats_grid(Db *db,
                                  Db *dbgrid,
                                  const char *oper,
                                  int ncol,
                                  int *cols,
                                  int radius)
{
  int icol, jcol, iptr, iptn, iptm, ndim, count;
  int *indg, *indg0, nxyz, error, i, iad, iech, nmed;
  double *coor, *medtab, ratio, value, mean, valdef;

  /* Initializations */

  error = 1;
  iptm = iptn = 0;
  indg = indg0 = nullptr;
  coor = medtab = nullptr;
  nxyz = dbgrid->getSampleNumber();
  ndim = dbgrid->getNDim();
  count = (int) pow(2. * radius + 1., (double) ndim);

  /* Preliminary check */

  if (!is_grid(dbgrid))
  {
    messerr("The Output Db must be a Grid");
    return (1);
  }

  /* Check that all variables are defined */

  for (icol = 0; icol < ncol; icol++)
  {
    jcol = cols[icol];
    if (!db->isColumnIndexValid(jcol))
    {
      messerr("Error: Variable %d is not defined", cols[icol]);
      goto label_end;
    }
  }

  /* Check the validity of the requested function */

  if (!st_oper_check(oper, 0, 1, 0, 1, 0)) goto label_end;

  /* Create and initialize the new attributes */

  valdef = 0.;
  if (!strcmp(oper, "mini")) valdef = 1.e30;
  if (!strcmp(oper, "maxi")) valdef = -1.e30;

  /* Create the attributes */

  if (!strcmp(oper, "mean") || !strcmp(oper, "var") || !strcmp(oper, "stdv"))
    iptn = dbgrid->addFieldsByConstant(1, 0.);
  if (!strcmp(oper, "var") || !strcmp(oper, "stdv"))
    iptm = dbgrid->addFieldsByConstant(1, 0.);

  /* Core allocation */

  indg = db_indg_alloc(dbgrid);
  if (indg == nullptr) goto label_end;
  indg0 = db_indg_alloc(dbgrid);
  if (indg0 == nullptr) goto label_end;
  coor = db_sample_alloc(db, ELoc::X);
  if (coor == nullptr) goto label_end;
  if (!strcmp(oper, "med"))
  {
    medtab = db_sample_alloc(db, ELoc::X);
    if (medtab == nullptr) goto label_end;
  }

  /* Loop on the variables */

  for (icol = 0; icol < ncol; icol++)
  {

    /* Create the output attribute */

    jcol = cols[icol];
    iptr = dbgrid->addFieldsByConstant(1, valdef);
    if (iptn > 0) db_attribute_init(dbgrid, 1, iptn, 0.);
    if (iptm > 0) db_attribute_init(dbgrid, 1, iptm, 0.);

    /* Loop on the samples */

    nmed = 0;
    for (iech = 0; iech < db->getSampleNumber(); iech++)
    {

      /* Read a sample */

      if (!db->isActive(iech)) continue;
      db_sample_load(db, ELoc::X, iech, coor);
      if (point_to_grid(dbgrid, coor, 0, indg0) < 0) continue;
      value = db->getArray(iech, jcol);
      if (FFFF(value)) continue;

      /* Loop on the neighboring cells */

      for (int ic = 0; ic < count; ic++)
      {
        st_get_neighboring_cell(ndim, radius, ic, indg0, indg);
        iad = db_index_grid_to_sample(dbgrid, indg);
        if (iad < 0 || iad >= nxyz) continue;

        if (!strcmp(oper, "num"))
        {
          dbgrid->updArray(iad, iptr, 0, 1.);
        }
        else if (!strcmp(oper, "mean"))
        {
          dbgrid->updArray(iad, iptn, 0, 1.);
          dbgrid->updArray(iad, iptr, 0, value);
        }
        else if (!strcmp(oper, "var") || !strcmp(oper, "stdv"))
        {
          dbgrid->updArray(iad, iptn, 0, 1.);
          dbgrid->updArray(iad, iptm, 0, value);
          dbgrid->updArray(iad, iptr, 0, value * value);
        }
        else if (!strcmp(oper, "mini"))
        {
          if (value < dbgrid->getArray(iad, iptr))
            dbgrid->setArray(iad, iptr, value);
        }
        else if (!strcmp(oper, "maxi"))
        {
          if (value > dbgrid->getArray(iad, iptr))
            dbgrid->setArray(iad, iptr, value);
        }
        else if (!strcmp(oper, "med"))
        {
          medtab[nmed++] = value;
        }
        else if (!strcmp(oper, "plus"))
        {
          if (value > 0.) dbgrid->updArray(iad, iptr, 0, 1.);
        }
        else if (!strcmp(oper, "moins"))
        {
          if (value < 0.) dbgrid->updArray(iad, iptr, 0, 1.);
        }
        else if (!strcmp(oper, "zero"))
        {
          if (value == 0.) dbgrid->updArray(iad, iptr, 0, 1.);
        }
        else
        {
          value = 0.;
        }
      }
    }

    /* Normation */

    if (!strcmp(oper, "mean"))
    {
      for (i = 0; i < nxyz; i++)
      {
        ratio = dbgrid->getArray(i, iptn);
        if (ratio <= 0.)
          dbgrid->setArray(i, iptr, TEST);
        else
          dbgrid->updArray(i, iptr, 3, ratio);
      }
    }
    else if (!strcmp(oper, "var") || !strcmp(oper, "stdv"))
    {
      for (i = 0; i < nxyz; i++)
      {
        ratio = dbgrid->getArray(i, iptn);
        if (ratio <= 0.)
          dbgrid->setArray(i, iptr, TEST);
        else
        {
          mean = dbgrid->getArray(i, iptm) / ratio;
          value = dbgrid->getArray(i, iptr) / ratio - mean * mean;
          if (value < 0) value = 0.;
          if (!strcmp(oper, "var"))
            dbgrid->setArray(i, iptr, value);
          else
            dbgrid->setArray(i, iptr, sqrt(value));
        }
      }
    }
    else if (!strcmp(oper, "med"))
    {
      for (i = 0; i < nxyz; i++)
      {
        value = (nmed > 0) ? medtab[nmed / 2] :
                             TEST;
        dbgrid->setArray(i, iptr, value);
      }
    }
    else if (!strcmp(oper, "mini"))
    {
      for (i = 0; i < nxyz; i++)
      {
        value = dbgrid->getArray(i, iptr);
        if (value == 1.e30) dbgrid->setArray(i, iptr, TEST);
      }
    }
    else if (!strcmp(oper, "maxi"))
    {
      for (i = 0; i < nxyz; i++)
      {
        value = dbgrid->getArray(i, iptr);
        if (value == -1.e30) dbgrid->setArray(i, iptr, TEST);
      }
    }
  }

  /* Set the error return flag */

  error = 0;

  label_end:

  /* Delete auxiliary attributes for local calculations */

  if ((!strcmp(oper, "mean") || !strcmp(oper, "var") || !strcmp(oper, "stdv")) && iptn
      > 0) dbgrid->deleteFieldByAttribute(iptn);
  if ((!strcmp(oper, "var") || !strcmp(oper, "stdv")) && iptm > 0)
    dbgrid->deleteFieldByAttribute(iptm);

  indg = db_indg_free(indg);
  indg0 = db_indg_free(indg0);
  coor = db_sample_free(coor);
  medtab = db_sample_free(medtab);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculates the statistics of points within cells of a grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid     Db for the grid
 ** \param[in]  db         Db for the points
 ** \param[in]  oper       Name of the operation
 ** \li                     "num"   : number of active samples
 ** \li                     "mean"  : Mean of the first variable
 ** \li                     "sum"   : Sum of the first variable
 ** \li                     "std"   : Standard deviation of the first variable
 ** \li                     "var"   : Variance of the first variable
 ** \li                     "mean2" : Mean of the second variable
 ** \li                     "sum2"  : Sum of the second variable
 ** \li                     "std2"  : Standard deviation of the second variable
 ** \li                     "var2"  : Variance of the second variable
 ** \li                     "cov"   : Covariance between two variables
 ** \li                     "corr"  : Correlation between two variables
 ** \li                     "mini"  : Minimum of the first variable
 ** \li                     "maxi"  : Maximum of the first variable
 ** \li                     "ore"   : Tonnage quantity
 ** \li                     "metal" : Metal quantity
 ** \param[in]  iatt       Rank of the first attribute
 ** \param[in]  jatt       Rank of the second attribute
 ** \param[in]  ncut       Number of cutoffs
 ** \param[in]  cuts       Array of cutoffs (optional)
 **
 ** \param[out]  tab       Output array (Dimension: 1 or ncut)
 **
 *****************************************************************************/
int stats_point_to_grid(Db *dbgrid,
                                        Db *db,
                                        const char *oper,
                                        int iatt,
                                        int jatt,
                                        int ncut,
                                        double *cuts,
                                        double *tab)
{
  int flag_s1, flag_v1, flag_s2, flag_v2, flag_v12, flag_mini, flag_maxi;
  int *indg, nxyz, error, iad, flag1, flag2, flag_denorm, flag_q, flag_t;
  double *coor, *nn, *s1, *s2, *v1, *v2, *v12, *mini, *maxi, *cutval, z1, z2,
      ratio;

  /* Initializations */

  error = 1;
  z1 = z2 = 0.;
  nn = s1 = s2 = v1 = v2 = v12 = mini = maxi = cutval = nullptr;
  indg = nullptr;
  coor = nullptr;
  nxyz = dbgrid->getSampleNumber();
  flag1 = flag2 = flag_denorm = flag_q = flag_t = 0;
  flag_s1 = flag_s2 = flag_v1 = flag_v2 = flag_v12 = flag_mini = flag_maxi = 0;

  /* Check the operator validity */

  if (!st_oper_check(oper, 1, 0, 1, 0, 1)) goto label_end;

  /* Set the relevant flags */

  if (!strcmp(oper, "num"))
    flag1 = 1;
  else if (!strcmp(oper, "mean"))
    flag1 = flag_s1 = 1;
  else if (!strcmp(oper, "sum"))
    flag1 = flag_s1 = flag_denorm = 1;
  else if (!strcmp(oper, "stdv"))
    flag1 = flag_s1 = flag_v1 = 1;
  else if (!strcmp(oper, "var"))
    flag1 = flag_s1 = flag_v1 = 1;
  else if (!strcmp(oper, "mean2"))
    flag2 = flag_s2 = 1;
  else if (!strcmp(oper, "sum2"))
    flag2 = flag_s2 = flag_denorm = 1;
  else if (!strcmp(oper, "stdv2"))
    flag2 = flag_s2 = flag_v2 = 1;
  else if (!strcmp(oper, "var2"))
    flag2 = flag_s2 = flag_v2 = 1;
  else if (!strcmp(oper, "cov"))
    flag2 = flag_s1 = flag_s2 = flag_v12 = 1;
  else if (!strcmp(oper, "corr"))
    flag2 = flag_s1 = flag_s2 = flag_v1 = flag_v2 = flag_v12 = 1;
  else if (!strcmp(oper, "mini"))
    flag1 = flag_mini = 1;
  else if (!strcmp(oper, "maxi"))
    flag1 = flag_maxi = 1;
  else if (!strcmp(oper, "ore"))
    flag1 = flag_t = 1;
  else if (!strcmp(oper, "metal"))
    flag1 = flag_q = 1;
  else
    return (1);

  /* Core allocation */

  coor = db_sample_alloc(db, ELoc::X);
  if (coor == nullptr) goto label_end;
  indg = db_indg_alloc(dbgrid);
  if (indg == nullptr) goto label_end;
  nn = (double*) mem_alloc(sizeof(double) * nxyz, 0);
  if (nn == nullptr) goto label_end;
  if (flag_s1)
  {
    s1 = (double*) mem_alloc(sizeof(double) * nxyz, 0);
    if (s1 == nullptr) goto label_end;
  }
  if (flag_s2)
  {
    s2 = (double*) mem_alloc(sizeof(double) * nxyz, 0);
    if (s2 == nullptr) goto label_end;
  }
  if (flag_v1)
  {
    v1 = (double*) mem_alloc(sizeof(double) * nxyz, 0);
    if (v1 == nullptr) goto label_end;
  }
  if (flag_v2)
  {
    v2 = (double*) mem_alloc(sizeof(double) * nxyz, 0);
    if (v2 == nullptr) goto label_end;
  }
  if (flag_v12)
  {
    v12 = (double*) mem_alloc(sizeof(double) * nxyz, 0);
    if (v12 == nullptr) goto label_end;
  }
  if (flag_mini)
  {
    mini = (double*) mem_alloc(sizeof(double) * nxyz, 0);
    if (mini == nullptr) goto label_end;
  }
  if (flag_maxi)
  {
    maxi = (double*) mem_alloc(sizeof(double) * nxyz, 0);
    if (maxi == nullptr) goto label_end;
  }
  if (flag_t || flag_q)
  {
    cutval = (double*) mem_alloc(sizeof(double) * nxyz * ncut, 0);
    if (cutval == nullptr) goto label_end;
  }

  /* Array initializations */

  for (int i = 0; i < nxyz; i++)
  {
    nn[i] = 0.;
    if (flag_s1) s1[i] = 0.;
    if (flag_s2) s2[i] = 0.;
    if (flag_v1) v1[i] = 0.;
    if (flag_v2) v2[i] = 0.;
    if (flag_v12) v12[i] = 0.;
    if (flag_mini) mini[i] = TEST;
    if (flag_maxi) maxi[i] = TEST;
    if (flag_t || flag_q) for (int icut = 0; icut < ncut; icut++)
      cutval[icut + i * ncut] = 0;
  }

  /* Loop on the samples */

  for (int iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Check the variable(s) */

    if (flag1)
    {
      z1 = db->getArray(iech, iatt);
      if (FFFF(z1)) continue;
    }
    if (flag2)
    {
      z2 = db->getArray(iech, jatt);
      if (FFFF(z2)) continue;
    }

    /* Check the location of the data in the grid */

    db_sample_load(db, ELoc::X, iech, coor);
    if (point_to_grid(dbgrid, coor, 0, indg) < 0) continue;

    iad = db_index_grid_to_sample(dbgrid, indg);
    if (iad < 0 || iad >= nxyz) continue;
    nn[iad]++;
    if (flag_s1) s1[iad] += z1;
    if (flag_s2) s2[iad] += z2;
    if (flag_v1) v1[iad] += z1 * z1;
    if (flag_v2) v2[iad] += z2 * z2;
    if (flag_v12) v12[iad] += z1 * z2;
    if (flag_mini)
    {
      if (FFFF(mini[iad]) || z1 < mini[iad]) mini[iad] = z1;
    }
    if (flag_maxi)
    {
      if (FFFF(maxi[iad]) || z1 > maxi[iad]) maxi[iad] = z1;
    }
    if (flag_t)
    {
      for (int icut = 0; icut < ncut; icut++)
        if (z1 >= cuts[icut]) cutval[icut + iech * ncut] += 1.;
    }
    if (flag_q)
    {
      for (int icut = 0; icut < ncut; icut++)
        if (z1 >= cuts[icut]) cutval[icut + iech * ncut] += z1;
    }
  }

  /* Normation */

  for (int i = 0; i < nxyz; i++)
  {
    ratio = nn[i];
    if (ratio <= 0)
    {
      if (flag_s1) s1[i] = TEST;
      if (flag_s2) s2[i] = TEST;
      if (flag_v1) v1[i] = TEST;
      if (flag_v2) v2[i] = TEST;
      if (flag_v12) v12[i] = TEST;
      if (flag_mini) mini[i] = TEST;
      if (flag_maxi) maxi[i] = TEST;
      if (flag_t || flag_q) for (int icut = 0; icut < ncut; icut++)
        cutval[icut + i * ncut] = TEST;
    }
    else
    {
      if (flag_s1) s1[i] /= ratio;
      if (flag_s2) s2[i] /= ratio;
      if (flag_v1)
      {
        v1[i] = v1[i] / ratio - s1[i] * s1[i];
        v1[i] = (v1[i] < 0.) ? 0. :
                               sqrt(v1[i]);
      }
      if (flag_v2)
      {
        v2[i] = v2[i] / ratio - s2[i] * s2[i];
        v2[i] = (v2[i] < 0.) ? 0. :
                               sqrt(v2[i]);
      }
      if (flag_v12)
      {
        v12[i] = v12[i] / ratio - s1[i] * s2[i];
      }
      if (flag_denorm)
      {
        if (flag_s1) s1[i] *= ratio;
        if (flag_s2) s2[i] *= ratio;
      }
      if (flag_t || flag_q)
      {
        for (int icut = 0; icut < ncut; icut++)
          cutval[icut + i * ncut] /= ratio;
      }
    }
  }

  /* Dispatch according to the type of result expected */

  for (int i = 0; i < nxyz; i++)
  {
    if (!strcmp(oper, "num"))
      tab[i] = nn[i];
    else if (!strcmp(oper, "mean"))
      tab[i] = s1[i];
    else if (!strcmp(oper, "sum"))
      tab[i] = s1[i];
    else if (!strcmp(oper, "stdv"))
      tab[i] = v1[i];
    else if (!strcmp(oper, "var"))
      tab[i] = v1[i] * v1[i];
    else if (!strcmp(oper, "mean2"))
      tab[i] = s2[i];
    else if (!strcmp(oper, "sum2"))
      tab[i] = s2[i];
    else if (!strcmp(oper, "stdv2"))
      tab[i] = v2[i];
    else if (!strcmp(oper, "var2"))
      tab[i] = v2[i] * v2[i];
    else if (!strcmp(oper, "cov"))
      tab[i] = v12[i];
    else if (!strcmp(oper, "corr"))
    {
      if (v1[i] > 0. && v2[i] > 0.) tab[i] = v12[i] / (v1[i] * v2[i]);
    }
    else if (!strcmp(oper, "mini"))
      tab[i] = mini[i];
    else if (!strcmp(oper, "maxi"))
      tab[i] = maxi[i];
    else if (!strcmp(oper, "ore"))
      for (int icut = 0; icut < ncut; icut++)
        tab[i + icut * nxyz] = cutval[icut + i * ncut];
    else if (!strcmp(oper, "metal"))
      for (int icut = 0; icut < ncut; icut++)
        tab[i + icut * nxyz] = cutval[icut + i * ncut];
    else
      goto label_end;
  }

  /* Set the error return flag */

  error = 0;

  label_end: indg = db_indg_free(indg);
  coor = db_sample_free(coor);
  nn = (double*) mem_free((char* ) nn);
  s1 = (double*) mem_free((char* ) s1);
  s2 = (double*) mem_free((char* ) s2);
  v1 = (double*) mem_free((char* ) v1);
  v2 = (double*) mem_free((char* ) v2);
  v12 = (double*) mem_free((char* ) v12);
  mini = (double*) mem_free((char* ) mini);
  maxi = (double*) mem_free((char* ) maxi);
  cutval = (double*) mem_free((char* ) cutval);
  return (error);
}

/****************************************************************************/
/*!
 **  Update the proportions
 **
 ** \param[in]  dbin       Db for the input grid
 ** \param[in]  indg       Array of grid indices
 ** \param[in]  nfacies    Number of facies
 **
 ** \param[out] prop       Array of proportions
 **
 *****************************************************************************/
static void st_update_prop(Db *dbin, int *indg, int nfacies, double *prop)
{
  int ifac;

  ifac = (int) dbin->getVariable(db_index_grid_to_sample(dbin, indg), 0);
  if (ifac < 1 || ifac > nfacies) return;
  prop[ifac - 1] += 1.;
}

/****************************************************************************/
/*!
 **  Update the transitions
 **
 ** \param[in]  dbin       Db for the input grid
 ** \param[in]  pos        Rank of the montee axis
 ** \param[in]  indg       Array of grid indices
 ** \param[in]  nfacies    Number of facies
 ** \param[in]  orient     Orientation
 **
 ** \param[out] trans      Array of transitions
 **
 *****************************************************************************/
static void st_update_trans(Db *dbin,
                            int pos,
                            int *indg,
                            int nfacies,
                            int orient,
                            double *trans)
{
  int ifac1, ifac2, jpos;

  jpos = indg[pos] + orient;
  if (jpos <= 0 || jpos >= dbin->getNX(pos)) return;
  ifac1 = (int) dbin->getVariable(db_index_grid_to_sample(dbin, indg), 0);
  indg[pos] += orient;
  ifac2 = (int) dbin->getVariable(db_index_grid_to_sample(dbin, indg), 0);
  indg[pos] -= orient;

  if (ifac1 < 1 || ifac1 > nfacies || ifac2 < 1 || ifac2 > nfacies) return;

  trans[(ifac1 - 1) * nfacies + (ifac2 - 1)] += 1.;
}

/****************************************************************************/
/*!
 **  Scale the proportions and store the proportions
 **
 ** \param[in]  dbout      Db for the output grid
 ** \param[in]  iptr       Writing pointer
 ** \param[in]  iech       Rank of the target sample
 ** \param[in]  nitem      Number of items
 ** \param[in]  tab        Array of cumulative statistics
 **
 *****************************************************************************/
static void st_scale_and_affect(Db *dbout,
                                int iptr,
                                int iech,
                                int nitem,
                                double *tab)
{
  int ifac;
  double value, total;

  total = 0.;
  for (ifac = 0; ifac < nitem; ifac++)
    total += tab[ifac];

  for (ifac = 0; ifac < nitem; ifac++)
  {
    if (total <= 0.)
      value = TEST;
    else
      value = tab[ifac] / total;
    dbout->setArray(iech, iptr + ifac, value);
  }
}

/****************************************************************************/
/*!
 **  Calculates the montee from a grid into a 1-D grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Db for the input grid
 ** \param[in]  dbout      Db for the output grid
 ** \param[in]  pos        Rank of the montee axis (starting from 0)
 ** \param[in]  nfacies    Number of facies
 ** \param[in]  radius     Radius of the neighborhood
 **
 *****************************************************************************/
int stats_proportion(Db *dbin,
                                     Db *dbout,
                                     int pos,
                                     int nfacies,
                                     int radius)
{
  double *prop;
  int *indg, ndim, ngrid, error, iptr, iech, ifac, aux1, aux2, bux1, bux2,
      ishift, i1, i2;

  /* Initializations */

  error = 1;
  indg = nullptr;
  prop = nullptr;

  /* Preliminary checks */

  ndim = dbin->getNDim();
  if (ndim != 2 && ndim != 3)
  {
    messerr("This function is limited to 2-D or 3-D input grids");
    goto label_end;
  }
  if (pos < 0 || pos >= ndim)
  {
    messerr("The rank of the 'montee' axis should lie between 1 and %d", ndim);
    goto label_end;
  }
  if (dbin->getNX(pos) != dbout->getNX(0) || dbin->getX0(pos) != dbout->getX0(0)
      || dbin->getDX(pos) != dbout->getDX(0))
  {
    messerr("The 1-D output grid does not match input grid");
    goto label_end;
  }
  if (!dbin->isVariableNumberComparedTo(1)) goto label_end;

  /* Core allocation */

  ngrid = dbin->getNX(pos);
  indg = db_indg_alloc(dbin);
  prop = (double*) mem_alloc(sizeof(double) * nfacies, 0);
  if (prop == nullptr) goto label_end;

  /* Create the new variables in the output file */

  iptr = dbout->addFieldsByConstant(nfacies, TEST);
  if (iptr < 0) goto label_end;

  /* Loop on the elements of the output grid */

  for (iech = 0; iech < ngrid; iech++)
  {
    for (ifac = 0; ifac < nfacies; ifac++)
      prop[ifac] = 0.;

    if (ndim == 2)
    {
      aux1 = (pos + 1) % ndim;
      for (ishift = -radius; ishift <= radius; ishift++)
      {
        indg[pos] = iech + ishift;
        if (indg[pos] < 0 || indg[pos] >= ngrid) continue;
        for (i1 = 0; i1 < dbin->getNX(aux1); i1++)
        {
          indg[aux1] = i1;
          st_update_prop(dbin, indg, nfacies, prop);
        }
        st_scale_and_affect(dbout, iptr, iech, nfacies, prop);
      }
    }
    else
    {
      bux1 = (pos + 1) % ndim;
      bux2 = (pos + 2) % ndim;
      aux1 = MIN(bux1, bux2);
      aux2 = MAX(bux1, bux2);
      for (ishift = -radius; ishift <= radius; ishift++)
      {
        indg[pos] = iech + ishift;
        if (indg[pos] < 0 || indg[pos] >= ngrid) continue;
        for (i1 = 0; i1 < dbin->getNX(aux1); i1++)
          for (i2 = 0; i2 < dbin->getNX(aux2); i2++)
          {
            indg[aux1] = i1;
            indg[aux2] = i2;
            st_update_prop(dbin, indg, nfacies, prop);
          }
        st_scale_and_affect(dbout, iptr, iech, nfacies, prop);
      }
    }
  }

  /* Set the error return flag */

  error = 0;

  label_end: prop = (double*) mem_free((char* ) prop);
  indg = db_indg_free(indg);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculates the transition from a grid into a 1-D grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       Db for the input grid
 ** \param[in]  dbout      Db for the output grid
 ** \param[in]  pos        Rank of the montee axis (starting from 0)
 ** \param[in]  nfacies    Number of facies
 ** \param[in]  radius     Radius of the neighborhood
 ** \param[in]  orient     Orientation (+1 or -1)
 **
 *****************************************************************************/
int stats_transition(Db *dbin,
                                     Db *dbout,
                                     int pos,
                                     int nfacies,
                                     int radius,
                                     int orient)
{
  double *trans;
  int *indg;
  int ndim, ngrid, error, iptr, iech, item, aux1, aux2, bux1, bux2, i1, i2,
      ishift, nitem;

  /* Initializations */

  error = 1;
  indg = nullptr;
  trans = nullptr;

  /* Preliminary checks */

  ndim = dbin->getNDim();
  if (ndim != 2 && ndim != 3)
  {
    messerr("This function is limited to 2-D or 3-D input grids");
    goto label_end;
  }
  if (pos < 0 || pos >= ndim)
  {
    messerr("The rank of the 'montee' axis should lie between 1 and %d", ndim);
    goto label_end;
  }
  if (dbin->getNX(pos) != dbout->getNX(0) || dbin->getX0(pos) != dbout->getX0(0)
      || dbin->getDX(pos) != dbout->getDX(0))
  {
    messerr("The 1-D output grid does not match input grid");
    goto label_end;
  }
  if (!dbin->isVariableNumberComparedTo(1)) goto label_end;

  /* Core allocation */

  ngrid = dbin->getNX(pos);
  indg = db_indg_alloc(dbin);
  nitem = nfacies * nfacies;
  trans = (double*) mem_alloc(sizeof(double) * nitem, 0);
  if (trans == nullptr) goto label_end;

  /* Create the new variables in the output file */

  iptr = dbout->addFieldsByConstant(nfacies * nfacies, TEST);
  if (iptr < 0) goto label_end;

  /* Loop on the elements of the output grid */

  for (iech = 0; iech < ngrid; iech++)
  {
    indg[pos] = iech;
    for (item = 0; item < nitem; item++)
      trans[item] = 0.;

    if (ndim == 2)
    {
      aux1 = (pos + 1) % ndim;
      for (ishift = -radius; ishift <= radius; ishift++)
      {
        indg[pos] = iech + ishift;
        if (indg[pos] < 0 || indg[pos] >= ngrid) continue;
        for (i1 = 0; i1 < dbin->getNX(aux1); i1++)
        {
          indg[aux1] = i1;
          st_update_trans(dbin, pos, indg, nfacies, orient, trans);
        }
        st_scale_and_affect(dbout, iptr, iech, nitem, trans);
      }
    }
    else
    {
      bux1 = (pos + 1) % ndim;
      bux2 = (pos + 2) % ndim;
      aux1 = MIN(bux1, bux2);
      aux2 = MAX(bux1, bux2);
      for (ishift = -radius; ishift <= radius; ishift++)
      {
        indg[pos] = iech + ishift;
        if (indg[pos] < 0 || indg[pos] >= ngrid) continue;
        for (i1 = 0; i1 < dbin->getNX(aux1); i1++)
          for (i2 = 0; i2 < dbin->getNX(aux2); i2++)
          {
            indg[aux1] = i1;
            indg[aux2] = i2;
            st_update_trans(dbin, pos, indg, nfacies, orient, trans);
          }
        st_scale_and_affect(dbout, iptr, iech, nitem, trans);
      }
    }
  }

  /* Set the error return flag */

  error = 0;

  label_end: trans = (double*) mem_free((char* ) trans);
  indg = db_indg_free(indg);
  return (error);
}

/****************************************************************************/
/*!
 **  Load the subgrid from the Input Db
 **
 ** \return Return the total probability of finding joins
 **
 ** \param[in]  verbose   Verbose flag
 ** \param[in]  flag_ffff 1 replace masked values by 0
 **                       0 replace masked values by FFFF
 ** \param[in]  iech0     Rank of the output grid cell
 ** \param[in]  nech0     Number of cells of the output grid
 ** \param[in]  ntot      Dimension of arrays 'numtab1' and 'valtab1'
 ** \param[in]  dbgrid    Db for the input grid
 ** \param[in]  ind0      Origin of the Output Db within the Input Db
 ** \param[in]  nxyz      Mesh of the Output Db (expressed in Input Db)
 ** \param[in]  ixyz      Indices of the starting node of the Output Db
 **
 ** \param[out] numtab1   Array containing the sample count
 ** \param[out] valtab1   Array containing the sample value
 **
 ** \remarks  The array ind0, ixyz and nxyz are dimensionned to ndim
 **
 *****************************************************************************/
static double st_extract_subgrid(int verbose,
                                 int flag_ffff,
                                 int iech0,
                                 int nech0,
                                 int ntot,
                                 Db *dbgrid,
                                 int *ind0,
                                 int *ixyz,
                                 int *nxyz,
                                 double *numtab1,
                                 double *valtab1)
{
  int ix, iy, iz, jx, jy, jz, ind, ecr, ndim;
  double proba, value;

  /* Initializations */

  ndim = dbgrid->getNDim();
  VectorInt iwork2(ndim);
  for (int i = 0; i < ntot; i++)
  {
    numtab1[i] = 0;
    valtab1[i] = 0.;
  }

  for (int idim = 0; idim < 3; idim++)
  {
    if (idim < ndim) continue;
    ixyz[idim] = 0;
    nxyz[idim] = 1;
    ind0[idim] = 0;
  }

  ecr = 0;
  proba = 0.;
  for (iz = 0; iz < nxyz[2]; iz++)
    for (iy = 0; iy < nxyz[1]; iy++)
      for (ix = 0; ix < nxyz[0]; ix++)
      {

        /* Get the address of a sample of the subgrid */

        jx = ind0[0] + ixyz[0] * nxyz[0] + ix;
        if (jx < 0 || jx > dbgrid->getNX(0)) continue;
        jy = ind0[1] + ixyz[1] * nxyz[1] + iy;
        if (jy < 0 || jy > dbgrid->getNX(1)) continue;
        jz = ind0[2] + ixyz[2] * nxyz[2] + iz;
        if (jz < 0 || jz > dbgrid->getNX(2)) continue;

        /* Get the node index within the Input Db */

        if (ndim >= 1) iwork2[0] = jx;
        if (ndim >= 2) iwork2[1] = jy;
        if (ndim >= 3) iwork2[2] = jz;
        ind = db_index_grid_to_sample(dbgrid, iwork2.data());
        numtab1[ecr] = 1.;
        value = dbgrid->isActive(ind) ? dbgrid->getVariable(ind, 0) :
                                        TEST;
        if (FFFF(value))
          valtab1[ecr] = (flag_ffff) ? 0 :
                                       TEST;
        else
        {
          valtab1[ecr] = value;
          proba += value;
        }
        ecr++;
      }

  /* Optional verbose option */

  if (verbose)
  {
    message("Output cell %3d/%3d = %d", iech0 + 1, nech0, nxyz[0]);
    for (int idim = 1; idim < ndim; idim++)
      message("x%d", nxyz[idim]);
    message(" cells of Input Grid (Proba=%lf)\n", proba);
  }
  return (proba);
}

/****************************************************************************/
/*!
 **  Divide by 2 in integer (upper rounded value)
 **
 ** \return  1 if the input value is larger than 2; 0 otherwise
 **
 ** \param[in]  nxyz      Dimensions of the subgrid
 ** \param[in]  orient    Rank of the target direction
 **
 *****************************************************************************/
static int st_divide_by_2(int *nxyz, int orient)
{
  int ival;

  ival = nxyz[orient];
  if (ival <= 1) return (0);

  ival = (int) floor((double) (ival + 1.) / 2);
  nxyz[orient] = ival;
  return (1);
}

/****************************************************************************/
/*!
 **  Perform the arithmetic mean
 **
 ** \param[in]  idim      Direction of calculation
 ** \param[in]  nxyz1     Dimensions of the initial subgrid
 ** \param[in]  nxyz2     Dimensions of the final subgrid
 **
 ** \param[out] numtab1   Array containing the sample count
 ** \param[out] numtab2   Array containing the sample count
 ** \param[out] valtab1   Array containing the sample value
 ** \param[out] valtab2   Array containing the sample value
 **
 *****************************************************************************/
static void st_mean_arith(int idim,
                          int *nxyz1,
                          int *nxyz2,
                          double *numtab1,
                          double *numtab2,
                          double *valtab1,
                          double *valtab2)
{
  int ix, iy, iz, ix1, ix2, iy1, iy2, iz1, iz2;

  for (iz = 0; iz < nxyz2[2]; iz++)
    for (iy = 0; iy < nxyz2[1]; iy++)
      for (ix = 0; ix < nxyz2[0]; ix++)
      {
        N2_TAB(ix,iy,iz) = V2_TAB(ix,iy,iz) = 0.;
        switch (idim)
        {
          case 0:
            ix1 = 2 * ix;
            ix2 = 2 * ix + 1;
            if (D1_TAB(ix1, iy, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix1, iy, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix1,iy,iz) * V1_TAB(ix1, iy, iz);
            }
            if (ix2 < nxyz1[0] && D1_TAB(ix2, iy, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix2, iy, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix2,iy,iz) * V1_TAB(ix2, iy, iz);
            }
            break;

          case 1:
            iy1 = 2 * iy;
            iy2 = 2 * iy + 1;
            if (D1_TAB(ix, iy1, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy1, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy1,iz) * V1_TAB(ix, iy1, iz);
            }
            if (iy2 < nxyz1[1] && D1_TAB(ix, iy2, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy2, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy2,iz) * V1_TAB(ix, iy2, iz);
            }
            break;

          case 2:
            iz1 = 2 * iz;
            iz2 = 2 * iz + 1;
            if (D1_TAB(ix, iy, iz1))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy, iz1);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy,iz1) * V1_TAB(ix, iy, iz1);
            }
            if (iz2 < nxyz1[2] && D1_TAB(ix, iy, iz2))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy, iz2);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy,iz2) * V1_TAB(ix, iy, iz2);
            }
            break;
        }
        V2_TAB(ix,iy,iz) =
            (N2_TAB(ix,iy,iz) > 0) ?
                                     V2_TAB(ix,iy,iz) / N2_TAB(ix, iy, iz) :
                                     TEST;
      }
  return;
}

/****************************************************************************/
/*!
 **  Perform the harmonic mean
 **
 ** \param[in]  idim      Direction of calculation
 ** \param[in]  nxyz1     Dimensions of the initial subgrid
 ** \param[in]  nxyz2     Dimensions of the final subgrid
 **
 ** \param[out] numtab1   Array containing the input  sample count
 ** \param[out] numtab2   Array containing the output sample count
 ** \param[out] valtab1   Array containing the input  sample value
 ** \param[out] valtab2   Array containing the output sample value
 **
 *****************************************************************************/
static void st_mean_harmo(int idim,
                          int *nxyz1,
                          int *nxyz2,
                          double *numtab1,
                          double *numtab2,
                          double *valtab1,
                          double *valtab2)
{
  int ix, iy, iz, ix1, ix2, iy1, iy2, iz1, iz2;

  for (iz = 0; iz < nxyz2[2]; iz++)
    for (iy = 0; iy < nxyz2[1]; iy++)
      for (ix = 0; ix < nxyz2[0]; ix++)
      {
        N2_TAB(ix,iy,iz) = V2_TAB(ix,iy,iz) = 0.;
        switch (idim)
        {
          case 0:
            ix1 = 2 * ix;
            ix2 = 2 * ix + 1;
            if (D1_TAB(ix1, iy, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix1, iy, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix1,iy,iz) / V1_TAB(ix1, iy, iz);
            }
            if (ix2 < nxyz1[0] && D1_TAB(ix2, iy, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix2, iy, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix2,iy,iz) / V1_TAB(ix2, iy, iz);
            }
            break;

          case 1:
            iy1 = 2 * iy;
            iy2 = 2 * iy + 1;
            if (D1_TAB(ix, iy1, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy1, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy1,iz) / V1_TAB(ix, iy1, iz);
            }
            if (iy2 < nxyz1[1] && D1_TAB(ix, iy2, iz))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy2, iz);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy2,iz) / V1_TAB(ix, iy2, iz);
            }
            break;

          case 2:
            iz1 = 2 * iz;
            iz2 = 2 * iz + 1;
            if (D1_TAB(ix, iy, iz1))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy, iz1);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy,iz1) / V1_TAB(ix, iy, iz1);
            }
            if (iz2 < nxyz1[2] && D1_TAB(ix, iy, iz2))
            {
              N2_TAB(ix,iy,iz) += N1_TAB(ix, iy, iz2);
              V2_TAB(ix,iy,iz) += N1_TAB(ix,iy,iz2) / V1_TAB(ix, iy, iz2);
            }
            break;
        }
        V2_TAB(ix,iy,iz) =
            (ABS(V2_TAB(ix,iy,iz)) > 1.e-10) ?
                                               N2_TAB(ix,iy,iz) / V2_TAB(ix, iy,
                                                                         iz) :
                                               TEST;
      }
  return;
}

/****************************************************************************/
/*!
 **  Update the dimensions and calculate the current number of cells
 **
 ** \return Number of cells
 **
 ** \param[in]  nxyz1     Dimensions of the initial subgrid
 ** \param[in]  numtab1   Array containing the sample count
 ** \param[in]  valtab1   Array containing the sample value
 **
 ** \param[out] nxyz2     Dimensions of the final subgrid
 ** \param[out] numtab2   Array containing the sample count
 ** \param[out] valtab2   Array containing the sample value
 **
 *****************************************************************************/
static int st_recopy(int *nxyz1,
                     double *numtab1,
                     double *valtab1,
                     int *nxyz2,
                     double *numtab2,
                     double *valtab2)
{
  int i, ncell;

  /* Update the dimension */

  ncell = 1;
  for (i = 0; i < 3; i++)
  {
    nxyz2[i] = nxyz1[i];
    ncell *= nxyz1[i];
  }

  /* Update the contents of the arrays */

  for (i = 0; i < ncell; i++)
  {
    numtab2[i] = numtab1[i];
    valtab2[i] = valtab1[i];
  }

  return (ncell);
}

/****************************************************************************/
/*!
 **  Print the generated grids (for counts and values)
 **
 ** \param[in]  subtitle  Subtitle
 ** \param[in]  nxyz      Dimensions of the grid
 **
 ** \param[out] numtab    Array containing the sample count
 ** \param[out] valtab    Array containing the sample value
 **
 ****************************************************************************/
static void st_print_grid(const char *subtitle,
                          int nxyz[3],
                          double *numtab,
                          double *valtab)
{
  char string[100];
  int iz, shift;

  /* Initializations */

  shift = nxyz[0] * nxyz[1];

  /* Loop on the third dimension */

  for (iz = 0; iz < nxyz[2]; iz++)
  {
    (void) gslSPrintf(string, "%s Values (iz=%d)\n", subtitle, iz + 1);
    message(string);
    print_matrix(NULL, 0, 0, nxyz[0], nxyz[1], NULL, &valtab[iz * shift]);
    (void) gslSPrintf(string, "%s Counts (iz=%d)\n", subtitle, iz + 1);
    message(string);
    print_matrix(NULL, 0, 0, nxyz[0], nxyz[1], NULL, &numtab[iz * shift]);
  }
  message("\n");
}

/****************************************************************************/
/*!
 **  Print statistics after basic upscaling operation
 **
 ** \param[in]  title     Title for the printout
 ** \param[in]  nxyz      Grid dimension
 ** \param[out] valtab    Array containing the sample value
 **
 *****************************************************************************/
static void st_print_upscale(const char *title, int *nxyz, double *valtab)
{
  double mini, maxi, value;
  int lec, ndef;

  mini = 1.e30;
  maxi = -1.e30;
  lec = ndef = 0;
  for (int iz = 0; iz < nxyz[2]; iz++)
    for (int iy = 0; iy < nxyz[1]; iy++)
      for (int ix = 0; ix < nxyz[0]; ix++)
      {
        value = valtab[lec++];
        if (FFFF(value)) continue;
        if (value < mini) mini = value;
        if (value > maxi) maxi = value;
        ndef++;
      }

  message("%11s ", title);
  message("(%3d %3d %3d) : ", nxyz[0], nxyz[1], nxyz[2]);
  message("%lf %lf (Def=%d)\n", mini, maxi, ndef);
}

/****************************************************************************/
/*!
 **  Upscale the variable defined on a subgrid into one value
 **
 ** \param[in]  orient    Upscaling orientation (0 to 2)
 ** \param[in]  nxyz      Dimensions of the subgrid
 ** \param[in]  flag_save If the array must be saved for keypair
 **
 ** \param[out] numtab0   Array containing the sample count
 ** \param[out] numtab1   Array containing the sample count
 ** \param[out] numtab2   Array containing the sample count
 ** \param[out] valtab0   Array containing the sample value
 ** \param[out] valtab1   Array containing the sample value
 ** \param[out] valtab2   Array containing the sample value
 ** \param[out] res1      First result (Harmonic first)
 ** \param[out] res2      First result (Arithmetic first)
 **
 *****************************************************************************/
static void st_upscale(int orient,
                       int *nxyz,
                       int flag_save,
                       double *numtab0,
                       double *numtab1,
                       double *numtab2,
                       double *valtab0,
                       double *valtab1,
                       double *valtab2,
                       double *res1,
                       double *res2)
{
  int idim, nxyz1[3], nxyz2[3], ncell, flag_debug;

  /* Initializations */

  flag_debug = OptDbg::query(EDbg::UPSCALE);

  /**************************************/
  /* Getting the minimum upscaled value */
  /**************************************/

  if (flag_debug)
  {
    mestitle(1, "Looking for the Minimum Upscaled Value");
    st_print_grid("Initial", nxyz, numtab0, valtab0);
  }
  ncell = st_recopy(nxyz, numtab0, valtab0, nxyz1, numtab1, valtab1);
  ncell = st_recopy(nxyz, numtab0, valtab0, nxyz2, numtab2, valtab2);
  if (flag_save) st_print_upscale("Initial", nxyz1, valtab1);

  while (ncell > 1)
  {

    /* Harmonic mean in the flow direction */

    if (st_divide_by_2(nxyz2, orient))
    {
      st_mean_harmo(orient, nxyz1, nxyz2, numtab1, numtab2, valtab1, valtab2);
      ncell = st_recopy(nxyz2, numtab2, valtab2, nxyz1, numtab1, valtab1);
      if (flag_save) st_print_upscale("Harmonic", nxyz1, valtab1);
      if (flag_debug) st_print_grid("Harmonic", nxyz1, numtab1, valtab1);
    }

    /* Arithmetic mean orthogonal to the flow direction */

    for (idim = 0; idim < 3; idim++)
    {
      if (idim == orient) continue;
      if (st_divide_by_2(nxyz2, idim))
      {
        st_mean_arith(idim, nxyz1, nxyz2, numtab1, numtab2, valtab1, valtab2);
        ncell = st_recopy(nxyz2, numtab2, valtab2, nxyz1, numtab1, valtab1);
        if (flag_save) st_print_upscale("Arithmetic", nxyz1, valtab1);
        if (flag_debug) st_print_grid("Arithmetic", nxyz1, numtab1, valtab1);
      }
    }
  }
  *res1 = valtab1[0];

  /**************************************/
  /* Getting the maximum upscaled value */
  /**************************************/

  if (flag_debug)
  {
    mestitle(1, "Looking for the Maximum Upscaled Value\n");
    st_print_grid("Initial", nxyz, numtab0, valtab0);
  }
  ncell = st_recopy(nxyz, numtab0, valtab0, nxyz1, numtab1, valtab1);
  ncell = st_recopy(nxyz, numtab0, valtab0, nxyz2, numtab2, valtab2);
  if (flag_save) st_print_upscale("Initial", nxyz1, valtab1);

  while (ncell > 1)
  {

    /* Arithmetic mean orthogonal to the flow direction */

    for (idim = 0; idim < 3; idim++)
    {
      if (idim == orient) continue;
      if (st_divide_by_2(nxyz2, idim))
      {
        st_mean_arith(idim, nxyz1, nxyz2, numtab1, numtab2, valtab1, valtab2);
        ncell = st_recopy(nxyz2, numtab2, valtab2, nxyz1, numtab1, valtab1);
        if (flag_save) st_print_upscale("Arithmetic", nxyz1, valtab1);
        if (flag_debug) st_print_grid("Arithmetic", nxyz1, numtab1, valtab1);
      }
    }

    /* Harmonic mean in the flow direction */

    if (st_divide_by_2(nxyz2, orient))
    {
      st_mean_harmo(orient, nxyz1, nxyz2, numtab1, numtab2, valtab1, valtab2);
      ncell = st_recopy(nxyz2, numtab2, valtab2, nxyz1, numtab1, valtab1);
      if (flag_save) st_print_upscale("Harmonic", nxyz1, valtab1);
      if (flag_debug) st_print_grid("Harmonic", nxyz1, numtab1, valtab1);
    }
  }
  *res2 = valtab1[0];

  /* Final result obtained using geometric mean */

  return;
}

/****************************************************************************/
/*!
 **  Check if the first grid is a subgrid of the second one
 **
 ** \return  1 if grids are multiple
 **
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  title      Title of the work (used only for verbose case)
 ** \param[in]  dbgrid1    Db for the input grid
 ** \param[in]  dbgrid2    Db for the output grid
 **
 ** \param[out] ind0       Starting address
 ** \param[out] nxyz       Number of subdivisions
 ** \param[out] ntot       Total number of grid nodes
 **
 *****************************************************************************/
static int st_is_subgrid(int verbose,
                         const char *title,
                         Db *dbgrid1,
                         Db *dbgrid2,
                         int *ind0,
                         int *nxyz,
                         int *ntot)
{
  double d;
  int ndim;

  /* Initializations */

  ndim = dbgrid1->getNDim();

  (*ntot) = 1;
  for (int idim = 0; idim < ndim; idim++)
  {
    /* Initializations */

    ind0[idim] = 0;
    nxyz[idim] = 1;

    /* Is origin of the output grid is located on a node of input grid */

    d = (dbgrid2->getX0(idim) - dbgrid1->getX0(idim)) / dbgrid1->getDX(idim);
    if (!isInteger(d))
    {
      messerr(
          "The origin of the Output Grid does not coincide with a node of the Input Grid");
      return (0);
    }
    ind0[idim] = (int) floor(d + 0.5);

    /* Are grid meshes multiple */

    d = dbgrid2->getDX(idim) / dbgrid1->getDX(idim);
    if (!isInteger(d))
    {
      messerr(
          "The grid cell of the Output Grid is not a multiple of the grid cell of the Input Grid");
      return (0);
    }
    nxyz[idim] = (int) floor(d + 0.5);
    (*ntot) *= nxyz[idim];
  }

  if (verbose)
  {
    mestitle(1, title);
    message("- Number of Cells =");
    for (int idim = 0; idim < ndim; idim++)
      message(" %d", nxyz[idim]);
    message("\n");
    message("- Index of Origin =");
    for (int idim = 0; idim < ndim; idim++)
      message(" %d", ind0[idim]);
    message("\n");
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Upscale one variable from a grid Db into another grid Db
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid1    Db for the input grid
 ** \param[in]  dbgrid2    Db for the output grid
 ** \param[in]  orient     Upscaling direction (0 to 2)
 ** \param[in]  verbose    Verbose flag
 **
 *****************************************************************************/
int db_upscale(Db *dbgrid1,
                               Db *dbgrid2,
                               int orient,
                               int verbose)
{
  double *valtab0, *valtab1, *valtab2, *numtab0, *numtab1, *numtab2;
  double result1, result2, result, probtot;
  int error, ndim, ind0[3], nxyz[3], ixyz[3], iech, iptr, ntot, ncol;
  int flag_save, iech_save;

  /* Initializations */

  valtab0 = valtab1 = valtab2 = numtab0 = numtab1 = numtab2 = nullptr;
  error = 1;
  flag_save = 0;
  iech_save = (int) get_keypone("Upscale.Converge.Block", 0);

  /* Preliminary checks */

  ndim = dbgrid1->getNDim();
  if (ndim < 1 || ndim > 3)
  {
    messerr("This function is limited to 2-D or 3-D input grids");
    goto label_end;
  }
  if (!dbgrid1->isVariableNumberComparedTo(1)) goto label_end;
  if (orient < 1 || orient > ndim)
  {
    messerr("Inconsistency between Orientation (%d) and Space dimension (%d)",
            orient, ndim);
    goto label_end;
  }
  orient--;

  /* Check that the output grid is a subgrid of the input grid */

  if (!st_is_subgrid(verbose, "Upscaling", dbgrid1, dbgrid2, ind0, nxyz, &ntot))
    goto label_end;

  /* Create the new variable in the output file */

  ncol = 3;
  iptr = dbgrid2->addFieldsByConstant(ncol, TEST);
  if (iptr < 0) goto label_end;
  dbgrid2->setLocatorsByAttribute(ncol, iptr, ELoc::Z);

  /* Core allocation */

  numtab0 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (numtab0 == nullptr) goto label_end;
  numtab1 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (numtab1 == nullptr) goto label_end;
  numtab2 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (numtab2 == nullptr) goto label_end;
  valtab0 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (valtab0 == nullptr) goto label_end;
  valtab1 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (valtab1 == nullptr) goto label_end;
  valtab2 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (valtab2 == nullptr) goto label_end;

  /* Loop on the cells of the Output Grid */

  for (iech = 0; iech < dbgrid2->getSampleNumber(); iech++)
  {
    result1 = result2 = result = TEST;
    OptDbg::setIndex(iech + 1);
    flag_save = (iech == iech_save - 1);
    if (dbgrid2->isActive(iech))
    {
      db_index_sample_to_grid(dbgrid2, iech, ixyz);

      /* Load the subgrid to be upscaled */

      probtot = st_extract_subgrid(verbose, 0, iech, dbgrid2->getSampleNumber(),
                                   ntot, dbgrid1, ind0, ixyz, nxyz, numtab0,
                                   valtab0);

      if (probtot > 0)
      {

        /* Upscale the corresponding subgrid of the Input Grid */

        st_upscale(orient, nxyz, flag_save, numtab0, numtab1, numtab2, valtab0,
                   valtab1, valtab2, &result1, &result2);
        result = sqrt(result1 * result2);
      }
      else
      {
        result = result1 = result2 = TEST;
      }
    }

    /* Store the result */

    dbgrid2->setVariable(iech, 0, result1);
    dbgrid2->setVariable(iech, 1, result2);
    dbgrid2->setVariable(iech, 2, result);
  }

  /* Set the error return code */

  error = 0;

  label_end: OptDbg::setIndex(0);
  numtab0 = (double*) mem_free((char* ) numtab0);
  numtab1 = (double*) mem_free((char* ) numtab1);
  numtab2 = (double*) mem_free((char* ) numtab2);
  valtab0 = (double*) mem_free((char* ) valtab0);
  valtab1 = (double*) mem_free((char* ) valtab1);
  valtab2 = (double*) mem_free((char* ) valtab2);
  return (error);
}

/****************************************************************************/
/*!
 **  Calculate the euclidean distance between current and initial seed locations
 **
 ** \return Calculated squared distance
 **
 ** \param[in]  orient    Diffusion orientation (0 or the space rank dimension)
 ** \param[in]  ndim      Space dimension
 ** \param[in]  locini    Initial seed positions
 ** \param[in]  loccur    Current seed positions
 **
 *****************************************************************************/
static double st_squared_distance(int orient,
                                  int ndim,
                                  int *locini,
                                  int *loccur)
{
  double delta, dist;

  dist = 0.;
  for (int idim = 0; idim < ndim; idim++)
  {
    if (orient != 0 && orient != (idim + 1)) continue;
    delta = locini[idim] - loccur[idim];
    dist += delta * delta;
  }
  return (dist);
}

/****************************************************************************/
/*!
 **  Converts from sample index into grid indices
 **
 ** \param[in]  ndim  Space dimension
 ** \param[in]  ntot  Total number of cells in subgrid
 ** \param[in]  nxyz  Dimensions of the subgrid
 ** \param[in]  iech  Rank of the sample
 **
 ** \param[out]  indg Grid indices
 **
 *****************************************************************************/
static void st_sample_to_grid(int ndim,
                              int ntot,
                              int *nxyz,
                              int iech,
                              int *indg)
{
  for (int idim = ndim - 1; idim >= 0; idim--)
  {
    ntot /= nxyz[idim];
    indg[idim] = iech / ntot;
    iech -= indg[idim] * ntot;
  }
}

/****************************************************************************/
/*!
 **  Converts from grid indices into sample index
 **
 ** \return The absolute index
 **
 ** \param[in]  ndim  Space dimension
 ** \param[in]  nxyz  Dimensions of the subgrid
 ** \param[in]  indg  Grid indices
 **
 *****************************************************************************/
static int st_grid_to_sample(int ndim, int *nxyz, int *indg)
{
  int idim, ival;

  ival = indg[ndim - 1];
  if (ival < 0 || ival >= nxyz[ndim - 1]) return (-1);
  for (idim = ndim - 2; idim >= 0; idim--)
  {
    if (indg[idim] < 0 || indg[idim] >= nxyz[idim]) return (-1);
    ival = ival * nxyz[idim] + indg[idim];
  }
  return (ival);
}

/****************************************************************************/
/*!
 **  Find the location of the cell (within the possible cells) which is the
 **  closest to the target cell provided as argument
 **
 ** \return Rank of the cell
 **
 ** \param[in]  ntot  Number of possibilities
 ** \param[in]  tab   Array containing the sample value
 ** \param[in]  cell  Cell location
 **
 *****************************************************************************/
static int st_fixed_position(int ntot, double *tab, int cell)
{
  int j;

  for (int i = 0; i < ntot; i++)
  {
    j = cell + i;
    if (j < ntot && tab[j] > 0) return (j);
    j = cell - i;
    if (j >= 0 && tab[j] > 0) return (j);
  }
  return (ntot - 1);
}

/****************************************************************************/
/*!
 **  Find the cell linked to the probability
 **
 ** \return Rank of the cell
 **
 ** \param[in]  ntot  Number of possibilities
 ** \param[in]  tab   Array containing the sample value
 ** \param[in]  proba Local probability
 **
 *****************************************************************************/
static int st_find_cell(int ntot, double *tab, double proba)
{
  double sum1, sum2;

  sum1 = sum2 = 0.;
  for (int i = 0; i < ntot; i++)
  {
    sum2 += tab[i];
    if (proba >= sum1 && proba < sum2) return (i);
    sum1 = sum2;
  }
  return (ntot - 1);
}

/****************************************************************************/
/*!
 **  Migrate the seed to one of the available neighboring cells
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in]  n_nbgh    Number of neighboring cells
 ** \param[in]  nxyz      Dimensions of the subgrid
 ** \param[in]  nbgh      Array giving the neighboring cell location
 ** \param[in]  valwrk    Working array (Dimension: n_nbgh)
 ** \param[in]  valtab0   Array containing the sample value
 ** \param[in]  locwrk    Working array for shifted seed positions
 ** \param[in,out] loccur    Current seed position
 **
 *****************************************************************************/
static void st_migrate_seed(int ndim,
                            int n_nbgh,
                            int *nxyz,
                            int *nbgh,
                            double *valwrk,
                            double *valtab0,
                            int *locwrk,
                            int *loccur)
{
  int iabs, ivois;
  double probtot, proba;

  /* Count the number of available neighboring cells */

  probtot = 0.;
  for (ivois = 0; ivois < n_nbgh; ivois++)
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      locwrk[idim] = loccur[idim] + NBGH(ivois, idim);
      locwrk[idim] = get_mirror_sample(nxyz[idim], locwrk[idim]);
    }
    iabs = st_grid_to_sample(ndim, nxyz, locwrk);
    valwrk[ivois] = valtab0[iabs];
    probtot += valwrk[ivois];
  }

  /* Draw a migration at random */

  if (probtot > 0)
  {
    proba = law_uniform(0., probtot);
    ivois = st_find_cell(n_nbgh, valwrk, proba);
    for (int idim = 0; idim < ndim; idim++)
      loccur[idim] += NBGH(ivois, idim);
  }
}

/****************************************************************************/
/*!
 **  Print the position (debug option)
 **
 ** \param[in]  ndim      Space dimension
 ** \param[in]  iseed     Rank of the trajectory (from 0)
 ** \param[in]  iter      Rank of the iteration (from 0)
 ** \param[in]  tab       Array containing the coordinates of the position
 **
 ** \remark  This very verbose option only functions if verbose is TRUE
 ** \remark  and the internal flag DEBUG is TRUE (this requires compiling)
 **
 *****************************************************************************/
static void st_print_position(int ndim, int iseed, int iter, int *tab)
{
  if (!DEBUG) return;
  message("Trajectory %d - Iteration %d:", iseed + 1, iter + 1);
  for (int idim = 0; idim < ndim; idim++)
    message(" %4d", tab[idim]);
  message("\n");
}

/****************************************************************************/
/*!
 **  Calculate the average squared distance as a function of the iteration
 **
 ** \param[in]  orient    Diffusion orientation (0 or the space rank dimension)
 ** \param[in]  ndim      Space dimension
 ** \param[in]  ntot      Total number of cells in the subgrid
 ** \param[in]  nseed     Number of seeds
 ** \param[in]  niter     Number of iterations
 ** \param[in]  n_nbgh    Number of neighboring cells
 ** \param[in]  flag_save If the array must be saved for keypair
 ** \param[in]  probtot   Total probability of joins
 ** \param[in]  nxyz      Dimensions of the subgrid
 ** \param[in]  nbgh      Array giving the neighboring cell location
 ** \param[in]  tabini    Array of initial seed positions
 **                       Dimension: nseed * ndim
 ** \param[in]  tabcur    Array of current seed positions
 **                       Dimension: nseed * ndim
 ** \param[in]  tabwrk    Array of shifted seed positions
 **                       Dimension: ndim
 ** \param[in]  valwrk    Working array
 **                       Dimension: n_nbgh
 ** \param[in]  valtab0   Array containing the sample value
 ** \param[in]  verbose   Verbose option
 **
 ** \param[out] cvdist2   Array of squared distances by iteration
 ** \param[out] trsave    Array of trajectories (saved only if defined)
 **
 ** \remarks If the starting position is fixed, it is specified using:
 **              set_keypair("Fixed_Position",...)
 **
 *****************************************************************************/
static void st_updiff(int orient,
                      int ndim,
                      int ntot,
                      int nseed,
                      int niter,
                      int n_nbgh,
                      int flag_save,
                      double probtot,
                      int *nxyz,
                      int *nbgh,
                      int *tabini,
                      int *tabcur,
                      int *tabwrk,
                      double *valwrk,
                      double *valtab0,
                      int verbose,
                      double *cvdist2,
                      double *trsave)
{
  double d2, dmoy, proba;
  int rank, fixed_position, flag_fixed;

  /* Check if a fixed starting position has been defined */

  fixed_position = (int) get_keypone("Fixed_Position", -1);
  flag_fixed = fixed_position >= 0;

  /* Draw initial seed locations */

  for (int iseed = 0; iseed < nseed; iseed++)
  {
    if (flag_fixed)
    {
      rank = st_fixed_position(ntot, valtab0, fixed_position);
    }
    else
    {
      proba = law_uniform(0., probtot);
      rank = st_find_cell(ntot, valtab0, proba);
    }
    st_sample_to_grid(ndim, ntot, nxyz, rank, &TABINI(iseed, 0));
    for (int idim = 0; idim < ndim; idim++)
      TABCUR(iseed,idim) = TABINI(iseed, idim);
    if (verbose) st_print_position(ndim, iseed, -1, &TABCUR(iseed, 0));
  }

  /* Loop on the iterations */

  for (int iter = 0; iter < niter; iter++)
  {
    dmoy = 0.;

    /* Loop on the seed points */

    for (int iseed = 0; iseed < nseed; iseed++)
    {

      /* Migrate the seed */

      st_migrate_seed(ndim, n_nbgh, nxyz, nbgh, valwrk, valtab0, tabwrk,
                      &TABCUR(iseed, 0));

      /* Optional printout */

      if (verbose) st_print_position(ndim, iseed, iter, &TABCUR(iseed, 0));

      /* Calculate the distance between current and initial positions */

      d2 = st_squared_distance(orient, ndim, &TABINI(iseed, 0),
                               &TABCUR(iseed, 0));

      /* Save the trajectory (optional) */

      if (flag_save && trsave != nullptr)
      {
        for (int idim = 0; idim < ndim; idim++)
          TRAJEC(iseed,iter,idim) = TABCUR(iseed, idim);
      }

      /* Update the mean distance */

      dmoy += d2;
    }
    cvdist2[iter] = dmoy / (double) nseed;
  }

  return;
}

/****************************************************************************/
/*!
 **  Update the quantities needed for calculating the linear regression
 **  amongst a set of 2-D points
 **
 ** \param[in]      x         X value
 ** \param[in]      y         Y value
 **
 ** \param[in,out]  count     Number of samples
 ** \param[in,out]  sum_x     Sum of the X variable
 ** \param[in,out]  sum_y     Sum of the Y variable
 ** \param[in,out]  sum_xx    Sum of the X*X variable
 ** \param[in,out]  sum_xy    Sum of the X*Y variable
 **
 *****************************************************************************/
static void st_update_regression(double x,
                                 double y,
                                 double *count,
                                 double *sum_x,
                                 double *sum_y,
                                 double *sum_xx,
                                 double *sum_xy)
{
  (*count) += 1.;
  (*sum_x) += x;
  (*sum_y) += y;
  (*sum_xx) += x * x;
  (*sum_xy) += x * y;
}

/****************************************************************************/
/*!
 **  Derive the diffusion factor from the serie of squared distance as
 **  a function of time
 **
 ** \return The diffusion coefficient
 **
 ** \param[in]  niter      Number of iterations
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  pmid       Percentage of niter to store convergence
 ** \param[in]  flag_save  If the array must be saved for keypair
 ** \param[in]  cvdist2    Array containing the squared distances
 **
 ** \param[out] cvsave     Array containing the storage (Dimension: 3 * niter)
 **
 *****************************************************************************/
static double st_get_diff_coeff(int niter,
                                int verbose,
                                double pmid,
                                int flag_save,
                                double *cvdist2,
                                double *cvsave)
{
  double slope, origin, slope_ref, origin_ref, sum_x, sum_y, sum_xy, sum_xx,
      count;
  double mx, my, var, cov;
  int iter, rank_mid;

  slope_ref = origin_ref = TEST;
  rank_mid = (int) (pmid * niter / 100.);
  count = sum_x = sum_y = sum_xx = sum_xy = origin = slope = 0.;

  for (int jter = 0; jter < niter; jter++)
  {
    iter = niter - jter - 1;

    /* Calculate the average slope */

    st_update_regression((double) (iter + 1), cvdist2[iter], &count, &sum_x,
                         &sum_y, &sum_xx, &sum_xy);
    if (count > 1)
    {
      mx = sum_x / count;
      my = sum_y / count;
      var = sum_xx / count - mx * mx;
      cov = sum_xy / count - mx * my;
      slope = cov / var;
      origin = my - slope * mx;
      if (iter == rank_mid)
      {
        slope_ref = slope;
        origin_ref = origin;
      }
    }

    /* Optional printout */

    if (verbose && !FFFF(slope) && !FFFF(origin))
    {
      message("  Rank=%5d Slope=%lf Origin=%lf (Count=%d)", iter + 1, slope,
              origin, (int) count);
      if (iter == rank_mid) message(" - Stored");
      message("\n");
    }

    /* Optional storage */

    if (flag_save)
    {
      cvsave[3 * iter + 0] = cvdist2[iter];
      cvsave[3 * iter + 1] = slope;
      cvsave[3 * iter + 2] = origin;
    }
  }

  /* Optional printout */

  if (verbose) message("- Slope=%lf Origin=%lf\n", slope_ref, origin_ref);

  /* Optional saving */

  if (flag_save) set_keypair("Diffusion.Converge", 1, niter, 3, cvsave);

  return (slope_ref);
}

/****************************************************************************/
/*!
 **  Calculate the diffusion factor from a grid Db into another grid Db
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid1    Db for the input grid
 ** \param[in]  dbgrid2    Db for the output grid
 ** \param[in]  orient     Diffusion orientation (0 or the space rank dimension)
 ** \param[in]  niter      Number of iterations
 ** \param[in]  nseed      Number of seeds
 ** \param[in]  seed       Seed for the random number generation
 ** \param[in]  verbose    Verbose Option
 **
 ** \remarks The user may specify the type of storage
 ** \remarks 0: average squared distance; 1: slope; 2: origin
 ** \remarks      get.keypair("Diffusion.Converge.Option",0)
 ** \remarks The neighborhood is Block(0) or Cross(1) according to
 ** \remarks      get.keypair("Diffusion.Converge.Morpho",0)
 **
 ** \remarks The rank of the convergence if measured from a Starting to an Ending
 ** \remarks percentages of the 'niter' iterations. The 'Mid' percentage gives
 ** \remarks the value assigned to the output block
 ** \remarks      get.keypair("Diffusion.Converge.PMid",70)
 **
 ** \remarks The user can specify the block of interest by using:
 ** \remarks      get.keypair("Diffusion.Converge.Block")
 ** \remarks For the target block, the array of squared distances as a function
 ** \remarks of the iteration is returned using the keypair mechanism:
 ** \remarks      set.keypair("Diffusion.Converge")
 ** \remarks For the target block, the trajectories are saved using keypair;
 ** \remarks      set.keypair("Diffusion.Trajectory.XX")
 **
 *****************************************************************************/
int db_diffusion(Db *dbgrid1,
                                 Db *dbgrid2,
                                 int orient,
                                 int niter,
                                 int nseed,
                                 int seed,
                                 int verbose)
{
  double *valtab0, *numtab0, *valwrk, *cvdist2, *cvsave, *trsave;
  double diff_coeff, pmid, probtot;
  int error, ndim, ind0[3], nxyz[3], ixyz[3], iech, nech, iptr, opt_center;
  int ntot, iech_save, flag_save, opt_morpho, flag_traj;
  int *tabini, *tabcur, *tabwrk, *numrank, n_nbgh;
  char name[40];
  VectorInt nbgh;

  /* Initializations */

  valtab0 = valwrk = numtab0 = cvdist2 = cvsave = trsave = nullptr;
  tabini = tabcur = tabwrk = numrank = nullptr;
  n_nbgh = flag_save = 0;
  error = 1;
  iech_save = (int) get_keypone("Diffusion.Converge.Block", 0);
  opt_morpho = (int) get_keypone("Diffusion.Converge.Morpho", 1);
  opt_center = (int) get_keypone("Diffusion.Converge.Center", 1);
  flag_traj = (int) get_keypone("Diffusion.Flag.Trajectory", 0);
  pmid = get_keypone("Diffusion.Converge.PMid", 70.);
  if (seed != 0) law_set_random_seed(seed);

  /* Preliminary checks */

  ndim = dbgrid1->getNDim();
  nech = dbgrid2->getSampleNumber();
  if (ndim < 1 || ndim > 3)
  {
    messerr("This function is limited to 2-D or 3-D input grids");
    goto label_end;
  }
  if (!(orient == 0 || (orient >= 1 && orient <= ndim)))
  {
    messerr("Argument 'orient' (%d) can be 0 or one of the space dimension",
            orient);
    goto label_end;
  }
  if (!dbgrid1->isVariableNumberComparedTo(1)) goto label_end;
  if (pmid < 5 || pmid > 95)
  {
    messerr("'PMid' must lie between 5% and 95%");
    goto label_end;
  }

  /* Check that the output grid is a subgrid of the input grid */

  if (!st_is_subgrid(verbose, "Diffusion Coefficient", dbgrid1, dbgrid2, ind0,
                     nxyz, &ntot)) goto label_end;

  /* Core allocation */

  tabini = (int*) mem_alloc(sizeof(int) * ndim * nseed, 0);
  if (tabini == nullptr) goto label_end;
  tabcur = (int*) mem_alloc(sizeof(int) * ndim * nseed, 0);
  if (tabcur == nullptr) goto label_end;
  tabwrk = (int*) mem_alloc(sizeof(int) * ndim, 0);
  if (tabwrk == nullptr) goto label_end;
  numrank = (int*) mem_alloc(sizeof(int) * ndim * nseed, 0);
  if (numrank == nullptr) goto label_end;
  numtab0 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (numtab0 == nullptr) goto label_end;
  valtab0 = (double*) mem_alloc(sizeof(double) * ntot, 0);
  if (valtab0 == nullptr) goto label_end;
  cvdist2 = (double*) mem_alloc(sizeof(double) * niter, 0);
  if (cvdist2 == nullptr) goto label_end;
  cvsave = (double*) mem_alloc(sizeof(double) * niter * 3, 0);
  if (cvsave == nullptr) goto label_end;
  if (flag_traj)
  {
    trsave = (double*) mem_alloc(sizeof(double) * niter * nseed * ndim, 0);
    if (trsave == nullptr) goto label_end;
  }

  /* Allocate the neighboring displacement array */

  nbgh = gridcell_neigh(ndim, opt_morpho, 1, opt_center, verbose, &n_nbgh);
  valwrk = (double*) mem_alloc(sizeof(double) * n_nbgh, 0);
  if (valwrk == nullptr) goto label_end;

  /* Create the new variable in the output file */

  iptr = dbgrid2->addFieldsByConstant(1, TEST);
  if (iptr < 0) goto label_end;

  /* Loop on the cells of the Output Grid */

  for (iech = 0; iech < nech; iech++)
  {
    diff_coeff = TEST;
    OptDbg::setIndex(iech + 1);
    flag_save = (iech == iech_save - 1);
    if (dbgrid2->isActive(iech))
    {
      db_index_sample_to_grid(dbgrid2, iech, ixyz);

      /* Load the subgrid to be upscaled */

      probtot = st_extract_subgrid(verbose, 1, iech, nech, ntot, dbgrid1, ind0,
                                   ixyz, nxyz, numtab0, valtab0);

      if (probtot > 0)
      {

        /* Upscale the diffusion */

        st_updiff(orient, ndim, ntot, nseed, niter, n_nbgh, flag_save, probtot,
                  nxyz, nbgh.data(), tabini, tabcur, tabwrk, valwrk, valtab0,
                  verbose, cvdist2, trsave);

        /* Derive the diffusion coefficient */

        diff_coeff = st_get_diff_coeff(niter, verbose, pmid, flag_save, cvdist2,
                                       cvsave);

        /* Save the trajectory (optional) */

        if (flag_save && trsave != nullptr)
        {
          for (int iseed = 0; iseed < nseed; iseed++)
          {
            (void) gslSPrintf(name, "Diffusion.Trajectory.%d", iseed + 1);
            for (int iter = 0; iter < niter; iter++)
              for (int idim = 0; idim < ndim; idim++)
                TRAJEC(iseed,iter,idim) = dbgrid2->getCoordinate(iech, idim) +
                TRAJEC(iseed,iter,idim) * dbgrid1->getDX(idim);
            set_keypair(name, 1, niter, ndim, &TRAJEC(iseed, 0, 0));
          }
        }
      }
    }

    /* Store the result */

    dbgrid2->setArray(iech, iptr, diff_coeff);
  }

  /* Set the error return code */

  error = 0;

  label_end: OptDbg::setIndex(0);
  tabini = (int*) mem_free((char* ) tabini);
  tabcur = (int*) mem_free((char* ) tabcur);
  tabwrk = (int*) mem_free((char* ) tabwrk);
  numrank = (int*) mem_free((char* ) numrank);
  valwrk = (double*) mem_free((char* ) valwrk);
  numtab0 = (double*) mem_free((char* ) numtab0);
  valtab0 = (double*) mem_free((char* ) valtab0);
  cvdist2 = (double*) mem_free((char* ) cvdist2);
  cvsave = (double*) mem_free((char* ) cvsave);
  trsave = (double*) mem_free((char* ) trsave);
  return (error);
}

/****************************************************************************/
/*!
 **  Constitute the name of the row
 **
 ** \param[in]  radix       Radix for the different variables (optional)
 ** \param[in]  ncol        Number of variables
 ** \param[in]  icol        Rank of the variable
 ** \param[in]  name        Variables name
 ** \param[in]  string      String array
 **
 *****************************************************************************/
static void st_get_rowname(const String &radix,
                           int ncol,
                           int icol,
                           const String &name,
                           char *string)
{
  if (!radix.empty())
    (void) gslSPrintf(string, "%s-%d", radix.c_str(), icol + 1);
  else if (!name.empty())
    (void) gslSPrintf(string, "%s", name.c_str());
  else if (ncol > 1)
    (void) gslSPrintf(string, "Variable-%d", icol + 1);
  else
    (void) gslSPrintf(string, "Variable");
}

/****************************************************************************/
/*!
 **  Print the multivariate statistics between different variables of a Db
 **
 ** \param[in]  db          Db structure
 ** \param[in]  iatts_arg   Ranks of the attributes (empty = all)
 ** \param[in]  flag_iso    Restrain statistics to isotopic samples
 ** \param[in]  flag_correl 1 if the correlations must be calculated
 ** \param[in]  title       Title for the printout (optional)
 ** \param[in]  radix       Radix for the different variables (optional)
 ** \param[in]  opers       Array of operators
 **
 *****************************************************************************/
void db_stats_print(const Db *db,
                    const VectorInt &iatts_arg,
                    const VectorString &opers,
                    int flag_iso,
                    int flag_correl,
                    const String &title,
                    const String &radix)
{
  double *data, *mean, *var, *mini, *maxi, *cov, *num;
  int iech, icol, jcol, numiso, ijcol, nundef, taille, noper, ncol;
  char string[50];

  /* Initializations */

  data = mean = mini = maxi = var = cov = num = nullptr;
  noper = static_cast<int>(opers.size());
  VectorInt iatts = iatts_arg;
  if (iatts.empty()) iatts = db->getAllAttributes();
  ncol = static_cast<int>(iatts.size());

  /* Preliminary checks */

  if (noper > 0)
  {
    for (int i = 0; i < noper; i++)
    {
      if (!st_oper_check(opers[i], 0, 0, 0, 0, 0)) goto label_end;
    }
  }
  if (flag_correl && ncol <= 1)
    messerr("Correlation matrix will not be printed for a single variable");

  /* Core allocation */

  data = (double*) mem_alloc(sizeof(double) * ncol, 0);
  if (data == nullptr) goto label_end;
  mean = (double*) mem_alloc(sizeof(double) * ncol, 0);
  if (mean == nullptr) goto label_end;
  mini = (double*) mem_alloc(sizeof(double) * ncol, 0);
  if (mini == nullptr) goto label_end;
  maxi = (double*) mem_alloc(sizeof(double) * ncol, 0);
  if (maxi == nullptr) goto label_end;
  var = (double*) mem_alloc(sizeof(double) * ncol, 0);
  if (var == nullptr) goto label_end;
  num = (double*) mem_alloc(sizeof(double) * ncol, 0);
  if (num == nullptr) goto label_end;
  if (flag_correl)
  {
    cov = (double*) mem_alloc(sizeof(double) * ncol * ncol, 0);
    if (cov == nullptr) goto label_end;
  }

  /* Initializations */

  numiso = 0;
  for (icol = ijcol = 0; icol < ncol; icol++)
  {
    mean[icol] = var[icol] = num[icol] = 0.;
    mini[icol] = 1.e30;
    maxi[icol] = -1.e30;
    if (flag_correl) for (jcol = 0; jcol < ncol; jcol++, ijcol++)
      cov[ijcol] = 0.;
  }

  /* Loop on the samples */

  for (iech = 0; iech < db->getSampleNumber(); iech++)
  {
    if (!db->isActive(iech)) continue;

    /* Look for isotopic sample */

    for (icol = nundef = 0; icol < ncol; icol++)
    {
      data[icol] = db->getArray(iech, iatts[icol]);
      if (FFFF(data[icol])) nundef++;
    }
    if (flag_iso && nundef > 0) continue;

    /* Calculate the 1-point statistics */

    for (icol = 0; icol < ncol; icol++)
    {
      if (FFFF(data[icol])) continue;
      num[icol] += 1.;
      mean[icol] += data[icol];
      var[icol] += data[icol] * data[icol];
      if (data[icol] < mini[icol]) mini[icol] = data[icol];
      if (data[icol] > maxi[icol]) maxi[icol] = data[icol];
    }

    if (nundef > 0) continue;
    numiso++;
    if (flag_correl) for (icol = ijcol = 0; icol < ncol; icol++)
      for (jcol = 0; jcol < ncol; jcol++, ijcol++)
        cov[ijcol] += data[icol] * data[jcol];
  }

  /* Normalization */

  for (icol = 0; icol < ncol; icol++)
  {
    if (num[icol] > 0)
    {
      mean[icol] /= num[icol];
      var[icol] /= num[icol];
      var[icol] -= mean[icol] * mean[icol];
      if (var[icol] <= 0) var[icol] = 0.;
    }
  }
  if (numiso > 0 && flag_correl) for (icol = ijcol = 0; icol < ncol; icol++)
    for (jcol = 0; jcol < ncol; jcol++, ijcol++)
    {
      cov[ijcol] /= numiso;
      cov[ijcol] -= mean[icol] * mean[jcol];
      cov[ijcol] /= sqrt(var[icol] * var[jcol]);
    }

  /************/
  /* Printout */
  /************/

  if (!title.empty()) mestitle(1, title.c_str());

  /* Calculate the maximum size of the variable */

  taille = 0;
  for (icol = 0; icol < ncol; icol++)
  {
    st_get_rowname(radix, ncol, icol, db_name_get_by_att(db, iatts[icol]),
                   string);
    taille = MAX(taille, (int ) strlen(string));
  }

  /* Print the header of the monovariate statistics */

  tab_print_rowname(" ", taille);
  if (st_oper_exists(opers, "num"))
    tab_prints(NULL, 1, EJustify::RIGHT, "Number");
  if (st_oper_exists(opers, "mini"))
    tab_prints(NULL, 1, EJustify::RIGHT, "Minimum");
  if (st_oper_exists(opers, "maxi"))
    tab_prints(NULL, 1, EJustify::RIGHT, "Maximum");
  if (st_oper_exists(opers, "mean"))
    tab_prints(NULL, 1, EJustify::RIGHT, "Mean");
  if (st_oper_exists(opers, "stdv"))
    tab_prints(NULL, 1, EJustify::RIGHT, "St. Dev.");
  if (st_oper_exists(opers, "var"))
    tab_prints(NULL, 1, EJustify::RIGHT, "Variance");
  message("\n");

  /* Print the monovariate statistics */

  for (icol = 0; icol < ncol; icol++)
  {
    st_get_rowname(radix, ncol, icol, db_name_get_by_att(db, iatts[icol]),
                   string);
    tab_print_rowname(string, taille);

    if (st_oper_exists(opers, "num"))
      tab_printi(NULL, 1, EJustify::RIGHT, (int) num[icol]);
    if (num[icol] > 0)
    {
      if (st_oper_exists(opers, "mini"))
        tab_printg(NULL, 1, EJustify::RIGHT, mini[icol]);
      if (st_oper_exists(opers, "maxi"))
        tab_printg(NULL, 1, EJustify::RIGHT, maxi[icol]);
      if (st_oper_exists(opers, "mean"))
        tab_printg(NULL, 1, EJustify::RIGHT, mean[icol]);
      if (st_oper_exists(opers, "stdv"))
        tab_printg(NULL, 1, EJustify::RIGHT, sqrt(var[icol]));
      if (st_oper_exists(opers, "var"))
        tab_printg(NULL, 1, EJustify::RIGHT, var[icol]);
    }
    else
    {
      if (st_oper_exists(opers, "mini"))
        tab_prints(NULL, 1, EJustify::RIGHT, "NA");
      if (st_oper_exists(opers, "maxi"))
        tab_prints(NULL, 1, EJustify::RIGHT, "NA");
      if (st_oper_exists(opers, "mean"))
        tab_prints(NULL, 1, EJustify::RIGHT, "NA");
      if (st_oper_exists(opers, "stdv"))
        tab_prints(NULL, 1, EJustify::RIGHT, "NA");
      if (st_oper_exists(opers, "var"))
        tab_prints(NULL, 1, EJustify::RIGHT, "NA");
    }
    message("\n");
  }
  message("\n");

  /* Print the correlation matrix  and count of isotopic samples */

  if (ncol > 1 && numiso > 0 && flag_correl)
  {
    message("Number of isotopic active samples = %d\n", numiso);
    print_matrix("Correlation matrix", 0, 1, ncol, ncol, NULL, cov);
    message("\n");
  }

  label_end: data = (double*) mem_free((char* ) data);
  mean = (double*) mem_free((char* ) mean);
  mini = (double*) mem_free((char* ) mini);
  maxi = (double*) mem_free((char* ) maxi);
  var = (double*) mem_free((char* ) var);
  cov = (double*) mem_free((char* ) cov);
  num = (double*) mem_free((char* ) num);
  return;
}

void db_stats_print(const Db *db,
                                    const VectorString &names,
                                    const VectorString &opers,
                                    int flag_iso,
                                    int flag_correl,
                                    const String &title,
                                    const String &radix)
{
  VectorInt iatts = db->getAttributes(names);
  if (iatts.size() <= 0) return;
  db_stats_print(db, iatts, opers, flag_iso, flag_correl, title, radix);
}

/****************************************************************************/
/*!
 **  Create residuals
 **
 ** \return  Error returned code
 **
 ** \param[in]  verbose Verbose flag
 ** \param[in]  nech    Number of samples
 ** \param[in]  tab     Array of sample values (Dimension: nech)
 ** \param[in]  ncut    Number of cutoffs
 ** \param[in]  zcut    Array of cutoff values (Dimension: ncut)
 **
 ** \param[out] nsorted   Number of sorted samples
 ** \param[out] mean      Average of the active data
 ** \param[out] residuals Array of residuals (Dimension: ncut * nech)
 ** \param[out] T         Array of for tonnage
 ** \param[out] Q         Array of for metal quantity
 **
 *****************************************************************************/
int stats_residuals(int verbose,
                                    int nech,
                                    double *tab,
                                    int ncut,
                                    double *zcut,
                                    int *nsorted,
                                    double *mean,
                                    double *residuals,
                                    double *T,
                                    double *Q)
{
  double value, moyenne;
  int iech, icut, jcut, nactive;

  /* Initializations */

  nactive = (*nsorted) = 0;
  moyenne = 0.;
  for (icut = 0; icut < ncut; icut++)
  {
    T[icut] = Q[icut] = 0.;
    for (iech = 0; iech < nech; iech++)
      RESIDUALS(icut,iech) = 0.;
  }

  /* Loop on the samples to calculate the indicators */

  for (iech = 0; iech < nech; iech++)
  {
    value = tab[iech];
    if (FFFF(value)) continue;
    moyenne += value;
    nactive++;

    /* Loop on the cutoffs */

    for (icut = 0; icut < ncut; icut++)
    {
      if (value < zcut[icut]) continue;
      RESIDUALS(icut,iech) = 1.;
      Q[icut] += value;
      T[icut] += 1.;
    }
  }
  if (nactive <= 0)
  {
    messerr("The calculation failed as there is no active sample");
    return (1);
  }

  /* Calculate the tonnage and meal quantity per class */

  moyenne /= (double) nactive;
  for (icut = 0; icut < ncut; icut++)
  {
    T[icut] /= (double) nactive;
    Q[icut] /= (double) nactive;
  }

  /* Calculate the residuals */

  for (iech = 0; iech < nech; iech++)
  {
    value = tab[iech];
    if (FFFF(value)) continue;

    /* Loop on the cutoffs */

    for (icut = ncut - 1; icut >= 0; icut--)
    {
      value = RESIDUALS(icut,iech) / T[icut];
      if (icut > 0)
      {
        jcut = icut - 1;
        value -= RESIDUALS(jcut,iech) / T[jcut];
      }
      else
      {
        value -= 1.;
      }
      RESIDUALS(icut,iech) = value;
    }
  }

  /* Verbose optional option */

  if (verbose)
  {
    mestitle(0, "Building residuals");
    message("Number of sorted samples = %d\n", nactive);
    for (icut = 0; icut < ncut; icut++)
      message("Cutoff %2d (above %lf) - Tonnage = %lf - Metal = %lf\n",
              icut + 1, zcut[icut], T[icut], Q[icut]);
  }

  (*nsorted) = nactive;
  (*mean) = moyenne;
  return (0);
}

