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
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/Utilities.hpp"
#include "geoslib_d.h"

/*! \cond */
#define INTRESX(ic, i)        (ctables->res[ic][i])
#define COVAL(ctables, iconf) (ctables->cmin + iconf * ctables->dc)
#define NELEM(ctables)        ((ctables->flag_cumul) ? ctables->ndisc + 1 : ctables->ndisc)
/*! \endcond */

namespace gstlrn
{
/****************************************************************************
 **
 ** FUNCTION: st_tableone_manage
 **
 ** PURPOSE: Management of one item of the CTables structure
 **
 ** RETURNS: For Initialization: 1 if actually performed; 0 otherwise
 ** RETURNS: For Freeing: 1 if actually performed; 0 otherwise
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  mode    : 1 for allocation; -1 for deallocation
 ** IN_ARGS:  rank    : rank within the CTables structure
 **
 ** OUT_ARGS: nb_used : Number of defined items for this discretization level
 **                     0 if the discretization level has not been used
 ** OUT_ARGS: nb_max  : Number of pixels in the Discretized covariance array
 **
 *****************************************************************************/
static void st_tableone_manage(CTables* ctables,
                               Id mode,
                               Id rank,
                               Id* nb_used,
                               Id* nb_max)
{
  Id nelem, size, number;

  // Initializations

  nelem    = NELEM(ctables);
  size     = nelem * nelem;
  *nb_used = 0;
  *nb_max  = size;

  // Dispatch

  if (mode > 0)
  {

    // Allocation

    if (ctables->res[rank].empty())
    {
      ctables->res[rank].resize(size);
      for (Id i = 0; i < size; i++)
        INTRESX(rank, i) = TEST;
      return;
    }
  }
  else
  {
    if (ctables->res[rank].empty())
    {
      number = 0;
      for (Id i = 0; i < size; i++)
        if (FFFF(INTRESX(rank, i))) number++;
      ctables->res[rank].clear();
      *nb_used = number;
      return;
    }
  }
}

/****************************************************************************
 **
 ** FUNCTION: ct_tableone_covrank
 **
 ** PURPOSE: Returns the rank of the target Ctable
 **
 ** RETURNS: Rank of the Target Ctable Rank
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  cova    : Covariance value
 **
 ** OUT_ARGS: cround  : Round discretized covariance value
 **
 *****************************************************************************/
Id ct_tableone_covrank(const CTables* ctables,
                       double cova,
                       double* cround)
{
  double dc, ecart;
  Id iconf, nconf;

  nconf            = ctables->nconf;
  ecart            = (cova - ctables->cmin);
  dc               = ctables->dc;
  auto placeholder = (0.5 + ecart / dc);
  iconf            = static_cast<Id>(placeholder);

  if (iconf < 0) iconf = 0;
  if (iconf >= nconf) iconf = nconf - 1;

  *cround = COVAL(ctables, iconf);

  return (iconf);
}

/****************************************************************************
 **
 ** FUNCTION: ct_INTRES2
 **
 ** PURPOSE: Get the value of the Discretized Table for given Configuration
 ** PURPOSE: If not available, the integral is calculated on the fly
 ** PURPOSE: Limitated to the 2-D case
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  iconf0  : Rank of the Target Ctable Rank
 ** IN_ARGS:  idisc0  : Index along first dimension
 ** IN_ARGS:  jdisc0  : Index along second dimension
 **
 *****************************************************************************/
double ct_INTRES2(CTables* ctables,
                  Id iconf0,
                  Id idisc0,
                  Id jdisc0)
{
  double lower[2], upper[2], error, value, cova;
  Id infin[2], inform, nelem, iad, nb_used, nb_max;
  static double abseps = 1.e-8;
  static double releps = 0.;
  static Id maxpts     = 25000;

  // Check if integral has already been defined

  if (ctables->res[iconf0].empty())
    st_tableone_manage(ctables, 1, iconf0, &nb_used, &nb_max);

  // Dispatch

  nelem = NELEM(ctables);
  iad   = idisc0 + nelem * jdisc0;
  cova  = COVAL(ctables, iconf0);

  // Check if the value has already been calculated

  value = INTRESX(iconf0, iad);
  if (!FFFF(value)) return (value);

  // First call to this value, calculate it

  if (!ctables->flag_cumul)
  {

    // Pixelated case

    lower[0] = ctables->v[idisc0];
    upper[0] = ctables->v[idisc0 + 1];
    infin[0] = 2;
    if (lower[0] == THRESH_INF) infin[0] = 0;
    if (upper[0] == THRESH_SUP) infin[0] = 1;

    lower[1] = ctables->v[jdisc0];
    upper[1] = ctables->v[jdisc0 + 1];
    infin[1] = 2;
    if (lower[1] == THRESH_INF) infin[1] = 0;
    if (upper[1] == THRESH_SUP) infin[1] = 1;

    mvndst(2, lower, upper, infin, &cova, maxpts, abseps, releps, &error,
           &value, &inform);
    if (inform) messageAbort("Error in function 'mvndst'");
  }
  else
  {

    // Cumulative case

    lower[0] = THRESH_INF;
    lower[1] = THRESH_INF;

    upper[0] = ctables->v[idisc0];
    infin[0] = 0;
    if (upper[0] == THRESH_SUP) infin[0] = 1;

    upper[1] = ctables->v[jdisc0];
    infin[1] = 0;
    if (upper[1] == THRESH_SUP) infin[1] = 1;

    mvndst(2, lower, upper, infin, &cova, maxpts, abseps, releps, &error,
           &value, &inform);
    if (inform) messageAbort("Error in function 'mvndst'");
  }

  // Store it in the table

  INTRESX(iconf0, iad) = value;
  return (value);
}

/****************************************************************************
 **
 ** FUNCTION: ct_INTRES3
 **
 ** PURPOSE: Get the value of the Discretized Table for given Configuration
 ** PURPOSE: If not available, the integral is calculated on the fly
 ** PURPOSE: Limitated to the 3-D case
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  iconf0  : Rank of the Target Ctable Rank
 ** IN_ARGS:  idisc0  : Index along first dimension
 ** IN_ARGS:  jdisc0  : Index along second dimension
 ** IN_ARGS:  kdisc0  : Index along third dimension
 **
 *****************************************************************************/
double ct_INTRES3(CTables* ctables,
                  Id iconf0,
                  Id idisc0,
                  Id jdisc0,
                  Id kdisc0)
{
  double lower[3], upper[3], error, value, cova;
  Id infin[3], inform, nelem, iad, nb_used, nb_max;
  static double abseps = 1.e-8;
  static double releps = 0.;
  static Id maxpts     = 25000;

  // Check if integral has already been defined

  if (ctables->res[iconf0].empty())
    st_tableone_manage(ctables, 1, iconf0, &nb_used, &nb_max);

  // Dispatch

  nelem = NELEM(ctables);
  iad   = idisc0 + nelem * (jdisc0 + nelem * kdisc0);
  cova  = COVAL(ctables, iconf0);

  // Check if the value has already been calculated

  value = INTRESX(iconf0, iad);
  if (!FFFF(value)) return (value);

  // First call to this value, calculate it

  if (!ctables->flag_cumul)
  {

    // Pixelated case

    lower[0] = ctables->v[idisc0];
    upper[0] = ctables->v[idisc0 + 1];
    infin[0] = 2;
    if (lower[0] == THRESH_INF) infin[0] = 0;
    if (upper[0] == THRESH_SUP) infin[0] = 1;

    lower[1] = ctables->v[jdisc0];
    upper[1] = ctables->v[jdisc0 + 1];
    infin[1] = 2;
    if (lower[1] == THRESH_INF) infin[1] = 0;
    if (upper[1] == THRESH_SUP) infin[1] = 1;

    lower[2] = ctables->v[kdisc0];
    upper[2] = ctables->v[kdisc0 + 1];
    infin[2] = 2;
    if (lower[2] == THRESH_INF) infin[2] = 0;
    if (upper[2] == THRESH_SUP) infin[2] = 1;

    mvndst(3, lower, upper, infin, &cova, maxpts, abseps, releps, &error,
           &value, &inform);
    if (inform) messageAbort("Error in function 'mvndst'");
  }
  else
  {

    // Cumulative case

    upper[0] = ctables->v[idisc0];
    infin[0] = 0;
    if (upper[0] == THRESH_SUP) infin[0] = 1;

    upper[1] = ctables->v[jdisc0];
    infin[1] = 0;
    if (upper[1] == THRESH_SUP) infin[1] = 1;

    upper[2] = ctables->v[kdisc0];
    infin[2] = 0;
    if (upper[2] == THRESH_SUP) infin[2] = 1;

    mvndst(3, lower, upper, infin, &cova, maxpts, abseps, releps, &error,
           &value, &inform);
    if (inform) messageAbort("Error in function 'mvndst'");
  }
  // Store it Id the table

  INTRESX(iconf0, iad) = value;
  return (value);
}

/****************************************************************************
 **
 ** FUNCTION:  ct_tables_print
 **
 ** PURPOSE: Print the CTables structure
 **
 ** IN_ARGS:  ctables    : Address to CTables
 ** IN_ARGS:  flag_print : Verbose option (0, 1 or 2)
 **
 *****************************************************************************/
void ct_tables_print(CTables* ctables, Id flag_print)
{
  Id ndisc, nconf, nelem;

  ndisc = ctables->ndisc;
  nconf = ctables->nconf;
  nelem = NELEM(ctables);

  mestitle(0, "Precalculation of Gaussian integral");
  message("Number of Covariance Discretizations steps  = %d\n", nconf);
  message("Lower Bound of Covariance Discretization    = %lf\n", ctables->cmin);
  message("Upper Bound of Covariance Discretization    = %lf\n", ctables->cmax);
  message("Covariance Discretization Interval          = %lf\n", ctables->dc);
  if (ctables->flag_cumul)
    message("Storing the integral from -Inf to the threshold per class\n");
  else
    message("Storing the integral per discretized class\n");
  message("Covariance is discretized between %lf and %lf\n", ctables->cmin,
          ctables->cmax);

  message("\n");
  message("Number of Probability Discretizations       = %d\n", ndisc);
  if (!ctables->v.empty())
    print_matrix("List of Gaussian Thresholds", 0, 1, 1, ctables->ndisc + 1,
                 NULL, ctables->v.data());

  if (flag_print > 0)
  {
    mestitle(2, "List of the configurations already calculated");

    for (Id iconf = 0; iconf < ctables->nconf; iconf++)
    {
      if (ctables->res[iconf].empty()) continue;

      if (flag_print > 0)
        message("- Configuration %d/%d (Cov=%lf)\n", iconf + 1, ctables->nconf,
                COVAL(ctables, iconf));

      if (flag_print == 2)
        print_matrix(NULL, 0, 1, nelem, nelem, NULL, ctables->res[iconf].data());
    }
    message("\n");
  }
}

/****************************************************************************
 **
 ** FUNCTION:  ct_tables_manage
 **
 ** PURPOSE: Management of the CTables structure
 **
 ** RETURNS: Pointer to the newly allocated CTables structure
 **
 ** IN_ARGS:  mode        : 1 for allocation; -1 for deallocation
 ** IN_ARGS:  verbose     : Verbose flag
 ** IN_ARGS:  flag_cumul  : 1 for storing gauss integral from -inf
 ** IN_ARGS:                0 for storing gauss integral per pixel
 ** IN_ARGS:  nconf       : Number of configurations
 ** IN_ARGS:  ndisc       : Number of Discretization steps
 ** IN_ARGS:  cmin        : Minimum value allowed for the correlation
 ** IN_ARGS:  cmax        : Maximum value allowed for the correlation
 ** IN_ARGS:  ctables_old : Address to CTables to be freed
 **
 *****************************************************************************/
CTables* ct_tables_manage(Id mode,
                          Id verbose,
                          Id flag_cumul,
                          Id nconf,
                          Id ndisc,
                          double cmin,
                          double cmax,
                          CTables* ctables_old)
{
  CTables* ctables;
  Id n_used, nb_used, nb_max;

  /* Dispatch */

  if (mode > 0)
  {
    // Allocation

    if (verbose)
      message("Allocating CTables (%dx%d) for %d possible configurations\n",
              ndisc, ndisc, nconf);
    ctables             = new CTables;
    ctables->flag_cumul = flag_cumul;
    ctables->nconf      = nconf;
    ctables->ndisc      = ndisc;
    ctables->cmin       = cmin;
    ctables->cmax       = cmax;
    ctables->dc         = (ctables->cmax - ctables->cmin) / static_cast<double>(nconf - 1);
    ctables->dp         = 1. / static_cast<double>(ndisc);

    ctables->res.resize(ctables->nconf);
    for (Id iconf = 0; iconf < ctables->nconf; iconf++)
      ctables->res[iconf].clear();

    // Define the array of thresholds

    ctables->v.resize(ndisc + 1);
    ctables->v[0]     = THRESH_INF;
    ctables->v[ndisc] = THRESH_SUP;
    for (Id idisc = 0; idisc < ndisc; idisc++)
      ctables->v[idisc] = law_invcdf_gaussian(static_cast<double>(idisc) * ctables->dp);
  }
  else
  {

    // De-allocation

    ctables = ctables_old;
    if (ctables == nullptr) return (ctables);
    ctables->v.clear();
    if (verbose)
      message("Freeing CTables from %d configuration(s)\n", ctables->nconf);

    n_used = 0;
    for (Id iconf = 0; iconf < ctables->nconf; iconf++)
    {
      st_tableone_manage(ctables, -1, iconf, &nb_used, &nb_max);
      if (nb_used > 0)
      {
        if (verbose)
          message("Configuration %3d - Number of items used: %d/%d\n",
                  iconf + 1, nb_used, nb_max);
        n_used++;
      }
    }
    if (verbose) message("Total of configurations actually used: %d\n", n_used);
    delete ctables;
    ctables = nullptr;
  }
  return (ctables);
}

/****************************************************************************
 **
 ** FUNCTION: st_tableone_getrank
 **
 ** PURPOSE: Get the starting and ending ranks of the integration
 ** PURPOSE: in the array of thresholds
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  low     : Lower bound
 ** IN_ARGS:  up      : Upper bound
 **
 ** OUT_ARGS: indmin  : index of the first discretized inside the interval
 ** OUT_ARGS: indmax  : index of the first discretized outside the interval
 **
 ** REMARKS:  The arguments (imin,imax) are returned so that we can easily
 ** REMARKS:  use them in a loop
 **
 *****************************************************************************/
static void st_tableone_getrank(const CTables* ctables,
                                double low,
                                double up,
                                Id* indmin,
                                Id* indmax)
{
  double v1;
  Id nelem;

  // Initializations

  nelem = NELEM(ctables);

  *indmin = *indmax = -1;

  // Dispatch

  if (!ctables->flag_cumul)
  {

    // Pixelated case

    for (Id idisc = 0; idisc < nelem - 1; idisc++)
    {
      v1 = (ctables->v[idisc] + ctables->v[idisc + 1]) / 2.;
      if (v1 < low) continue;
      if (*indmin < 0) *indmin = idisc;
      if (v1 > up)
      {
        *indmax = idisc;
        break;
      }
    }
    if (*indmax < 0) *indmax = nelem - 1;
  }
  else
  {

    // Cumulative case

    for (Id idisc = 0; idisc < nelem; idisc++)
    {
      v1 = ctables->v[idisc];
      if (v1 < low) continue;
      if (*indmin < 0) *indmin = MAX(0, idisc - 1);
      if (v1 >= up)
      {
        *indmax = MIN(nelem - 1, idisc);
        break;
      }
    }
    if (*indmax < 0) *indmax = nelem;
  }
}

/****************************************************************************
 **
 ** FUNCTION: ct_tableone_calculate
 **
 ** PURPOSE: Evaluate the integral over a rectangle from the integral
 ** PURPOSE: values calculated at pixel size
 **
 ** RETURNS: Integral value
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  iconf0  : Rank of the relevant internal structure
 ** IN_ARGS:  lows    : Array of lower values
 ** IN_ARGS:  ups     : Array of upper values
 **
 ** REMARKS: The code has been designed for boundaries of facies areas
 ** REMAKRS: parallel to the main axes
 **
 *****************************************************************************/
double ct_tableone_calculate(CTables* ctables,
                             Id iconf0,
                             double* lows,
                             double* ups)
{
  double result;
  Id i1min, i1max, i2min, i2max;

  // Initializations

  result = 0;

  // Dispatch

  if (!ctables->flag_cumul)
  {

    // Pixelated case

    st_tableone_getrank(ctables, lows[0], ups[0], &i1min, &i1max);
    st_tableone_getrank(ctables, lows[1], ups[1], &i2min, &i2max);
    for (Id idisc = i1min; idisc < i1max; idisc++)
      for (Id jdisc = i2min; jdisc < i2max; jdisc++)
        result += ct_INTRES2(ctables, iconf0, idisc, jdisc);
  }
  else
  {

    // Cumulative case

    st_tableone_getrank(ctables, lows[0], ups[0], &i1min, &i1max);
    st_tableone_getrank(ctables, lows[1], ups[1], &i2min, &i2max);
    result = (ct_INTRES2(ctables, iconf0, i1max, i2max) - ct_INTRES2(ctables, iconf0, i1min, i2max) - ct_INTRES2(ctables, iconf0, i1max, i2min) + ct_INTRES2(ctables, iconf0, i1min, i2min));
  }
  return (result);
}

/****************************************************************************
 **
 ** FUNCTION: ct_tableone_calculate_by_rank
 **
 ** PURPOSE: Evaluate the integral over a rectangle from the integral
 ** PURPOSE: values calculated at pixel size
 **
 ** RETURNS: Integral value
 **
 ** IN_ARGS:  ctables : Pointer to the CTables structure
 ** IN_ARGS:  iconf0  : Rank of the relevant internal structure
 ** IN_ARGS:  rklows  : Array of indices of lower values
 ** IN_ARGS:  rkups   : Array of indices upper values
 **
 ** REMARKS: The arguments 'rklows' and ;rkups' are double although
 ** REMARKS: they designate ranks.
 **
 *****************************************************************************/
double ct_tableone_calculate_by_rank(CTables* ctables,
                                     Id iconf0,
                                     double* rklows,
                                     double* rkups)
{
  double result;

  // Initializations

  result = 0;

  // Dispatch

  if (!ctables->flag_cumul)
  {

    // Pixelated case

    for (Id idisc = static_cast<Id>(rklows[0]); idisc < static_cast<Id>(rkups[0]); idisc++)
      for (Id jdisc = static_cast<Id>(rklows[1]); jdisc < static_cast<Id>(rkups[1]); jdisc++)
        result += ct_INTRES2(ctables, iconf0, idisc, jdisc);
  }
  else
  {

    // Cumulative case

    result = (ct_INTRES2(ctables, iconf0, static_cast<Id>(rkups[0]), static_cast<Id>(rkups[1])) -
              ct_INTRES2(ctables, iconf0, static_cast<Id>(rklows[0]), static_cast<Id>(rkups[1])) -
              ct_INTRES2(ctables, iconf0, static_cast<Id>(rkups[0]), static_cast<Id>(rklows[1])) +
              ct_INTRES2(ctables, iconf0, static_cast<Id>(rklows[0]), static_cast<Id>(rklows[1])));
  }
  return (result);
}

/****************************************************************************
 **
 ** FUNCTION: ct_tableone_getrank_from_proba
 **
 ** PURPOSE: Get the starting and ending ranks of the integration
 ** PURPOSE: in the array of thresholds (starting from proportions)
 **
 ** RETURNS: Index for the bound in the discretization table
 **
 ** IN_ARGS:  ctables  : Pointer to the CTables structure
 ** IN_ARGS:  gaussian : Gaussian bound
 **
 *****************************************************************************/
Id ct_tableone_getrank_from_proba(CTables* ctables,
                                  double gaussian)
{
  double dp, proba, vmin, vmax;
  Id iad, nelem;

  // Initializations

  nelem = NELEM(ctables);
  dp    = ctables->dp;

  proba = law_cdf_gaussian(gaussian);

  iad  = static_cast<Id>(proba / dp);
  vmin = dp * iad;
  vmax = dp * (iad + 1);
  if (vmax - proba < proba - vmin) iad++;

  iad = MIN(nelem, iad);

  return (iad);
}
} // namespace gstlrn